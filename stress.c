#include "stress.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAX_LINE 256

struct Joint {
    int id;
    double x, y;
    int is_support;
    char releases[8]; // e.g. "MX", "MY"
};

struct Member {
    int id;
    int j1_id, j2_id;
    double az, iy, ix;
};

struct JointLoad {
    int jid;
    char axis; // 'Z', 'X', 'Y'
    double magnitude;
};

struct MemberLoad {
    int mid;
    char axis;
    char type; // 'U' for Uniform, 'C' for Conc
    double vals[4];
};

struct StressModel {
    Joint* joints;
    int num_joints;
    int max_joints;

    Member* members;
    int num_members;
    int max_members;

    JointLoad* joint_loads;
    int num_joint_loads;
    int max_joint_loads;

    MemberLoad* member_loads;
    int num_member_loads;
    int max_member_loads;

    double E, G;
    char type[32];

    double* displacements; // Array of 3*num_joints
    double* reactions;     // Array of 3*num_joints
};

StressModel* stress_create(void) {
    StressModel* model = calloc(1, sizeof(StressModel));
    if (!model) return NULL;
    model->max_joints = 100;
    model->joints = calloc(model->max_joints, sizeof(Joint));
    model->max_members = 200;
    model->members = calloc(model->max_members, sizeof(Member));
    model->max_joint_loads = 100;
    model->joint_loads = calloc(model->max_joint_loads, sizeof(JointLoad));
    model->max_member_loads = 200;
    model->member_loads = calloc(model->max_member_loads, sizeof(MemberLoad));
    if (!model->joints || !model->members || !model->joint_loads || !model->member_loads) {
        stress_free(model);
        return NULL;
    }
    return model;
}

static int ensure_capacity(void** arr, int* cap, int count, size_t elem_size) {
    if (count < *cap) return 0;
    int new_cap = (*cap > 0) ? (*cap * 2) : 16;
    if (new_cap <= count) new_cap = count + 1;
    void* tmp = realloc(*arr, (size_t)new_cap * elem_size);
    if (!tmp) return -1;
    *arr = tmp;
    *cap = new_cap;
    return 0;
}

static void trim(char* s) {
    char* p = s;
    size_t l = strlen(p);
    while (l > 0 && isspace((unsigned char)p[l - 1])) p[--l] = 0;
    while (*p && isspace((unsigned char)*p)) {
        p++;
        l--;
    }
    memmove(s, p, l + 1);
}

int stress_load_dat(StressModel* model, const char* filename) {
    if (!model || !filename) return -1;
    FILE* f = fopen(filename, "r");
    if (!f) return -1;

    char line[MAX_LINE];
    char section[64] = "";

    while (fgets(line, MAX_LINE, f)) {
        trim(line);
        if (line[0] == 0 || line[0] == '*') continue;

        if (strncmp(line, "TYPE", 4) == 0 && (line[4] == '\0' || isspace((unsigned char)line[4]))) {
            char t1[16] = "", t2[16] = "";
            int n_tok = sscanf(line + 4, "%15s %15s", t1, t2);
            if (n_tok == 2) {
                snprintf(model->type, sizeof(model->type), "%s %s", t1, t2);
            } else if (n_tok == 1) {
                snprintf(model->type, sizeof(model->type), "%s", t1);
            }
            continue;
        }
        if (strncmp(line, "CONSTANTS", 9) == 0 && (line[9] == '\0' || isspace((unsigned char)line[9]))) {
            double nu = 0.2;
            char* e_ptr = strstr(line, "E ");
            if (e_ptr) sscanf(e_ptr + 2, "%lf", &model->E);
            char* p_ptr = strstr(line, "POISSON");
            if (p_ptr) sscanf(p_ptr + 7, "%lf", &nu);
            char* g_ptr = strstr(line, "G ");
            if (g_ptr && sscanf(g_ptr + 2, "%lf", &model->G) == 1) {
                /* explicit G */
            } else {
                model->G = model->E / (2.0 * (1.0 + nu));
            }
            continue;
        }

        // Section markers
        if (strstr(line, "JOINT COOR")) { strcpy(section, "COOR"); continue; }
        if (strstr(line, "MEMB INCI")) { strcpy(section, "INCI"); continue; }
        if (strstr(line, "MEMBER PROPERTIES")) { strcpy(section, "PROP"); continue; }
        if (strstr(line, "JOINT RELEASES")) { strcpy(section, "RELE"); continue; }
        if (strstr(line, "JOINT LOADS")) { strcpy(section, "JLOAD"); continue; }
        if (strstr(line, "MEMBER LOADS")) { strcpy(section, "MLOAD"); continue; }
        if (strstr(line, "SOLVE")) break;

        // Skip non-numeric lines in sections (like header commands)
        if (!isdigit((unsigned char)line[0])) continue;

        if (strcmp(section, "COOR") == 0) {
            Joint j = {0};
            if (sscanf(line, "%d %lf %lf", &j.id, &j.x, &j.y) >= 3) {
                if (strstr(line, "SUPPORT")) j.is_support = 1;
                if (ensure_capacity((void**)&model->joints, &model->max_joints,
                                    model->num_joints, sizeof(Joint)) != 0) {
                    fclose(f);
                    return -1;
                }
                model->joints[model->num_joints++] = j;
            }
        } else if (strcmp(section, "INCI") == 0) {
            Member m = {0};
            if (sscanf(line, "%d %d %d", &m.id, &m.j1_id, &m.j2_id) == 3) {
                if (ensure_capacity((void**)&model->members, &model->max_members,
                                    model->num_members, sizeof(Member)) != 0) {
                    fclose(f);
                    return -1;
                }
                model->members[model->num_members++] = m;
            }
        } else if (strcmp(section, "PROP") == 0) {
            int start = 0, end = 0;
            double az=0, iy=0, ix=0;
            if (strstr(line, "THRU")) {
                if (sscanf(line, "%d THRU %d", &start, &end) != 2) continue;
            } else {
                if (sscanf(line, "%d", &start) != 1) continue;
                end = start;
            }
            char* p;
            if ((p = strstr(line, "AZ"))) sscanf(p+2, "%lf", &az);
            if ((p = strstr(line, "IY"))) sscanf(p+2, "%lf", &iy);
            if ((p = strstr(line, "IX"))) sscanf(p+2, "%lf", &ix);
            for (int i=0; i < model->num_members; i++) {
                if (model->members[i].id >= start && model->members[i].id <= end) {
                    model->members[i].az = az;
                    model->members[i].iy = iy;
                    model->members[i].ix = ix;
                }
            }
        } else if (strcmp(section, "RELE") == 0) {
            int jid = 0;
            if (sscanf(line, "%d", &jid) != 1) continue;
            for (int i=0; i < model->num_joints; i++) {
                if (model->joints[i].id == jid) {
                    char* rel = model->joints[i].releases;
                    size_t cap = sizeof(model->joints[i].releases);
                    if (strstr(line, "MOMENT X") && !strstr(rel, "MX") && strlen(rel) + 2 < cap) {
                        strcat(rel, "MX");
                    }
                    if (strstr(line, "MOMENT Y") && !strstr(rel, "MY") && strlen(rel) + 2 < cap) {
                        strcat(rel, "MY");
                    }
                }
            }
        } else if (strcmp(section, "JLOAD") == 0) {
            JointLoad jl = {0};
            char axis_str[8] = "";
            if (sscanf(line, "%d FORCE %7s %lf", &jl.jid, axis_str, &jl.magnitude) == 3 ||
                sscanf(line, "%d MOMENT %7s %lf", &jl.jid, axis_str, &jl.magnitude) == 3) {
                jl.axis = (char)toupper((unsigned char)axis_str[0]);
                if (ensure_capacity((void**)&model->joint_loads, &model->max_joint_loads,
                                    model->num_joint_loads, sizeof(JointLoad)) != 0) {
                    fclose(f);
                    return -1;
                }
                model->joint_loads[model->num_joint_loads++] = jl;
            }
        } else if (strcmp(section, "MLOAD") == 0) {
            int start = 0, end = 0;
            if (strstr(line, "THRU")) {
                if (sscanf(line, "%d THRU %d", &start, &end) != 2) continue;
            } else {
                if (sscanf(line, "%d", &start) != 1) continue;
                end = start;
            }
            MemberLoad ml = {0};
            char axis_str[16] = "", type_str[32] = "";
            char* f_ptr = strstr(line, "FORCE");
            if (f_ptr) {
                int consumed = 0;
                if (sscanf(f_ptr, "FORCE %15s %31s%n", axis_str, type_str, &consumed) < 2) {
                    continue;
                }
                ml.axis = axis_str[0];
                ml.type = (strstr(type_str, "UNIFORM")) ? 'U' : 'C';
                
                // Parse up to 4 values after the type string
                char* val_ptr = f_ptr + consumed;
                int v_idx = 0;
                char* token = strtok(val_ptr, " \t");
                while (token && v_idx < 4) {
                    if (isdigit((unsigned char)token[0]) || token[0] == '-' || token[0] == '.') {
                        ml.vals[v_idx++] = strtod(token, NULL);
                    }
                    token = strtok(NULL, " \t");
                }

                for (int m_id = start; m_id <= end; m_id++) {
                    ml.mid = m_id;
                    if (ensure_capacity((void**)&model->member_loads, &model->max_member_loads,
                                        model->num_member_loads, sizeof(MemberLoad)) != 0) {
                        fclose(f);
                        return -1;
                    }
                    model->member_loads[model->num_member_loads++] = ml;
                }
            }
        }
    }
    fclose(f);
    return 0;
}

static int solve_system(double* K, double* F, size_t n, double* U) {
    // Gaussian Elimination with Partial Pivoting
    for (size_t i = 0; i < n; i++) {
        // Pivot
        size_t pivot = i;
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(K[j * n + i]) > fabs(K[pivot * n + i])) pivot = j;
        }
        if (fabs(K[pivot * n + i]) <= 1e-12 || isnan(K[pivot * n + i])) {
            return -1;
        }
        // Swap rows
        if (pivot != i) {
            for (size_t k = i; k < n; k++) {
                double temp = K[i * n + k];
                K[i * n + k] = K[pivot * n + k];
                K[pivot * n + k] = temp;
            }
            double temp_f = F[i];
            F[i] = F[pivot];
            F[pivot] = temp_f;
        }

        // Eliminate
        for (size_t j = i + 1; j < n; j++) {
            double factor = K[j * n + i] / K[i * n + i];
            F[j] -= factor * F[i];
            for (size_t k = i; k < n; k++) {
                K[j * n + k] -= factor * K[i * n + k];
            }
        }
    }
    // Back substitution
    for (size_t i = n; i > 0; i--) {
        size_t idx = i - 1;
        double sum = 0;
        for (size_t j = idx + 1; j < n; j++) {
            sum += K[idx * n + j] * U[j];
        }
        if (fabs(K[idx * n + idx]) <= 1e-12 || isnan(K[idx * n + idx])) {
            return -1;
        }
        U[idx] = (F[idx] - sum) / K[idx * n + idx];
    }
    return 0;
}

static int find_joint_index(const StressModel* model, int jid) {
    for (int j = 0; j < model->num_joints; j++) {
        if (model->joints[j].id == jid) return j;
    }
    return -1;
}

int stress_solve(StressModel* model) {
    if (!model || model->num_joints <= 0) return -1;

    free(model->displacements);
    model->displacements = NULL;
    free(model->reactions);
    model->reactions = NULL;

    int n = model->num_joints;
    size_t dofs = (size_t)n * 3;
    double* K = calloc(dofs * dofs, sizeof(double));
    double* F = calloc(dofs, sizeof(double));
    if (!K || !F) {
        free(K);
        free(F);
        return -1;
    }

    // Assembly
    for (int i = 0; i < model->num_members; i++) {
        Member m = model->members[i];
        int idx1 = find_joint_index(model, m.j1_id);
        int idx2 = find_joint_index(model, m.j2_id);
        if (idx1 < 0 || idx2 < 0) {
            free(K); free(F);
            return -1;
        }
        const Joint* j1 = &model->joints[idx1];
        const Joint* j2 = &model->joints[idx2];
        double dx = j2->x - j1->x;
        double dy = j2->y - j1->y;
        double L = sqrt(dx*dx + dy*dy);
        if (L <= 1e-12) {
            free(K); free(F);
            return -1;
        }
        double c = dx / L;
        double s = dy / L;

        double E = model->E, G = model->G, Iy = m.iy, J = m.ix;
        double L2 = L*L, L3 = L*L*L;
        double k_bend[4][4] = {
            {12*E*Iy/L3, 6*E*Iy/L2, -12*E*Iy/L3, 6*E*Iy/L2},
            {6*E*Iy/L2, 4*E*Iy/L, -6*E*Iy/L2, 2*E*Iy/L},
            {-12*E*Iy/L3, -6*E*Iy/L2, 12*E*Iy/L3, -6*E*Iy/L2},
            {6*E*Iy/L2, 2*E*Iy/L, -6*E*Iy/L2, 4*E*Iy/L}
        };
        double kt = G * J / L;

        double kl[6][6] = {0};
        int b_idx[] = {0, 2, 3, 5};
        for(int r=0; r<4; r++) for(int col=0; col<4; col++) kl[b_idx[r]][b_idx[col]] = k_bend[r][col];
        kl[1][1] = kt; kl[1][4] = -kt; kl[4][1] = -kt; kl[4][4] = kt;

        double T[6][6] = {0};
        double Tn[3][3] = {{1,0,0},{0,c,s},{0,-s,c}};
        for(int r=0; r<3; r++) for(int col=0; col<3; col++) { T[r][col] = Tn[r][col]; T[r+3][col+3] = Tn[r][col]; }

        // kg = T' * kl * T
        double temp[6][6] = {0};
        for(int r=0; r<6; r++) for(int col=0; col<6; col++) for(int k=0; k<6; k++) temp[r][col] += kl[r][k] * T[k][col];
        double kg[6][6] = {0};
        for(int r=0; r<6; r++) for(int col=0; col<6; col++) for(int k=0; k<6; k++) kg[r][col] += T[k][r] * temp[k][col];

        int dof_idx[] = {idx1*3, idx1*3+1, idx1*3+2, idx2*3, idx2*3+1, idx2*3+2};
        for(int r=0; r<6; r++) for(int col=0; col<6; col++) K[dof_idx[r]*dofs + dof_idx[col]] += kg[r][col];
    }

    // Loads
    for (int i = 0; i < model->num_joint_loads; i++) {
        int idx = find_joint_index(model, model->joint_loads[i].jid);
        if (idx < 0) {
            free(K); free(F);
            return -1;
        }
        char axis = model->joint_loads[i].axis;
        if (axis == 'Z') F[idx*3] += model->joint_loads[i].magnitude;
        else if (axis == 'X') F[idx*3+1] += model->joint_loads[i].magnitude;
        else if (axis == 'Y') F[idx*3+2] += model->joint_loads[i].magnitude;
    }
    for (int i = 0; i < model->num_member_loads; i++) {
        MemberLoad ml = model->member_loads[i];
        Member* m = NULL;
        for(int j=0; j<model->num_members; j++) if(model->members[j].id == ml.mid) m = &model->members[j];
        if (!m) continue;
        int idx1 = find_joint_index(model, m->j1_id);
        int idx2 = find_joint_index(model, m->j2_id);
        if (idx1 < 0 || idx2 < 0) {
            free(K); free(F);
            return -1;
        }
        const Joint* j1 = &model->joints[idx1];
        const Joint* j2 = &model->joints[idx2];
        double dx = j2->x - j1->x;
        double dy = j2->y - j1->y;
        double L = sqrt(dx*dx + dy*dy);
        if (L <= 1e-12) {
            free(K); free(F);
            return -1;
        }
        double c = dx / L;
        double s = dy / L;
        if (ml.type == 'U') {
            double w = ml.vals[0];
            double fef_l[6] = {w*L/2, 0, w*L*L/12, w*L/2, 0, -w*L*L/12};
            double Tn[3][3] = {{1,0,0},{0,c,s},{0,-s,c}};
            for (int k=0; k<3; k++) {
                for (int d=0; d<3; d++) {
                    F[idx1*3+d] += Tn[k][d] * fef_l[k];
                    F[idx2*3+d] += Tn[k][d] * fef_l[k+3];
                }
            }
        } else if (ml.type == 'C') {
            double p = ml.vals[0];
            if (ml.vals[1] < 0.0 || ml.vals[1] > L) {
                free(K); free(F);
                return -1;
            }
            double a = (ml.vals[1] > 0.0) ? ml.vals[1] : L/2;
            double b = L - a;
            double fef_w1 = p * b*b * (3*a + b) / (L*L*L);
            double fef_m1 = p * a * b*b / (L*L);
            double fef_w2 = p * a*a * (a + 3*b) / (L*L*L);
            double fef_m2 = -p * a*a * b / (L*L);
            double fef_l[6] = {fef_w1, 0, fef_m1, fef_w2, 0, fef_m2};
            double Tn[3][3] = {{1,0,0},{0,c,s},{0,-s,c}};
            for (int k=0; k<3; k++) {
                for (int d=0; d<3; d++) {
                    F[idx1*3+d] += Tn[k][d] * fef_l[k];
                    F[idx2*3+d] += Tn[k][d] * fef_l[k+3];
                }
            }
        }
    }

    // Boundary Conditions & Solution
    int* free_dofs = malloc(dofs * sizeof(int));
    if (!free_dofs) {
        free(K); free(F);
        return -1;
    }
    size_t num_free = 0;
    for (int i = 0; i < n; i++) {
        if (!model->joints[i].is_support) {
            free_dofs[num_free++] = i*3;
            free_dofs[num_free++] = i*3+1;
            free_dofs[num_free++] = i*3+2;
        } else {
            if (strstr(model->joints[i].releases, "MX")) free_dofs[num_free++] = i*3+1;
            if (strstr(model->joints[i].releases, "MY")) free_dofs[num_free++] = i*3+2;
        }
    }

    if (num_free == 0) {
        free(K); free(F); free(free_dofs);
        return -1;
    }

    double* K_sub = calloc(num_free * num_free, sizeof(double));
    double* F_sub = calloc(num_free, sizeof(double));
    double* U_sub = calloc(num_free, sizeof(double));
    if (!K_sub || !F_sub || !U_sub) {
        free(K); free(F); free(free_dofs); free(K_sub); free(F_sub); free(U_sub);
        return -1;
    }
    for(size_t r=0; r<num_free; r++) {
        F_sub[r] = F[free_dofs[r]];
        for(size_t col=0; col<num_free; col++) K_sub[r*num_free + col] = K[free_dofs[r]*dofs + free_dofs[col]];
    }

    if (solve_system(K_sub, F_sub, num_free, U_sub) != 0) {
        free(K); free(F); free(free_dofs); free(K_sub); free(F_sub); free(U_sub);
        return -1;
    }

    model->displacements = calloc(dofs, sizeof(double));
    model->reactions = calloc(dofs, sizeof(double));
    if (!model->displacements || !model->reactions) {
        free(model->displacements); model->displacements = NULL;
        free(model->reactions); model->reactions = NULL;
        free(K); free(F); free(free_dofs); free(K_sub); free(F_sub); free(U_sub);
        return -1;
    }
    for(size_t i=0; i<num_free; i++) model->displacements[free_dofs[i]] = U_sub[i];

    for(size_t r=0; r<dofs; r++) {
        double r_sum = 0;
        for(size_t col=0; col<dofs; col++) r_sum += K[r*dofs+col] * model->displacements[col];
        model->reactions[r] = r_sum - F[r];
    }

    free(K); free(F); free(free_dofs); free(K_sub); free(F_sub); free(U_sub);
    return 0;
}

void stress_free(StressModel* model) {
    if (!model) return;
    free(model->joints);
    free(model->members);
    free(model->joint_loads);
    free(model->member_loads);
    free(model->displacements);
    free(model->reactions);
    free(model);
}

/* Joint getters */
int stress_joint_get_id(const Joint* joint) { return joint ? joint->id : 0; }
double stress_joint_get_x(const Joint* joint) { return joint ? joint->x : 0.0; }
double stress_joint_get_y(const Joint* joint) { return joint ? joint->y : 0.0; }
int stress_joint_get_is_support(const Joint* joint) { return joint ? joint->is_support : 0; }
const char* stress_joint_get_releases(const Joint* joint) { return joint ? joint->releases : NULL; }

/* Member getters */
int stress_member_get_id(const Member* member) { return member ? member->id : 0; }
int stress_member_get_j1_id(const Member* member) { return member ? member->j1_id : 0; }
int stress_member_get_j2_id(const Member* member) { return member ? member->j2_id : 0; }
double stress_member_get_az(const Member* member) { return member ? member->az : 0.0; }
double stress_member_get_iy(const Member* member) { return member ? member->iy : 0.0; }
double stress_member_get_ix(const Member* member) { return member ? member->ix : 0.0; }

/* JointLoad getters */
int stress_joint_load_get_jid(const JointLoad* load) { return load ? load->jid : 0; }
char stress_joint_load_get_axis(const JointLoad* load) { return load ? load->axis : '\0'; }
double stress_joint_load_get_magnitude(const JointLoad* load) { return load ? load->magnitude : 0.0; }

/* MemberLoad getters */
int stress_member_load_get_mid(const MemberLoad* load) { return load ? load->mid : 0; }
char stress_member_load_get_axis(const MemberLoad* load) { return load ? load->axis : '\0'; }
char stress_member_load_get_type(const MemberLoad* load) { return load ? load->type : '\0'; }
const double* stress_member_load_get_vals(const MemberLoad* load) { return load ? load->vals : NULL; }

/* StressModel getters */
const Joint* stress_model_get_joint(const StressModel* model, int index) {
    if (!model || index < 0 || index >= model->num_joints) return NULL;
    return &model->joints[index];
}
int stress_model_get_num_joints(const StressModel* model) { return model ? model->num_joints : 0; }
int stress_model_get_max_joints(const StressModel* model) { return model ? model->max_joints : 0; }

const Member* stress_model_get_member(const StressModel* model, int index) {
    if (!model || index < 0 || index >= model->num_members) return NULL;
    return &model->members[index];
}
int stress_model_get_num_members(const StressModel* model) { return model ? model->num_members : 0; }
int stress_model_get_max_members(const StressModel* model) { return model ? model->max_members : 0; }

const JointLoad* stress_model_get_joint_load(const StressModel* model, int index) {
    if (!model || index < 0 || index >= model->num_joint_loads) return NULL;
    return &model->joint_loads[index];
}
int stress_model_get_num_joint_loads(const StressModel* model) { return model ? model->num_joint_loads : 0; }
int stress_model_get_max_joint_loads(const StressModel* model) { return model ? model->max_joint_loads : 0; }

const MemberLoad* stress_model_get_member_load(const StressModel* model, int index) {
    if (!model || index < 0 || index >= model->num_member_loads) return NULL;
    return &model->member_loads[index];
}
int stress_model_get_num_member_loads(const StressModel* model) { return model ? model->num_member_loads : 0; }
int stress_model_get_max_member_loads(const StressModel* model) { return model ? model->max_member_loads : 0; }

double stress_model_get_E(const StressModel* model) { return model ? model->E : 0.0; }
double stress_model_get_G(const StressModel* model) { return model ? model->G : 0.0; }
const char* stress_model_get_type(const StressModel* model) { return model ? model->type : NULL; }

const double* stress_model_get_displacements(const StressModel* model) {
    return model ? model->displacements : NULL;
}
const double* stress_model_get_reactions(const StressModel* model) {
    return model ? model->reactions : NULL;
}

/* Generic array mutator helpers */
static int insert_item(void** arr, int* count, int* cap, int index,
                       const void* elem, size_t elem_size) {
    if (index < 0 || index > *count) return -1;
    if (ensure_capacity(arr, cap, *count, elem_size) != 0) return -1;
    char* base = (char*)(*arr);
    if (index < *count) {
        memmove(base + (size_t)(index + 1) * elem_size,
                base + (size_t)index * elem_size,
                (size_t)(*count - index) * elem_size);
    }
    memcpy(base + (size_t)index * elem_size, elem, elem_size);
    (*count)++;
    return 0;
}

static int delete_item(void* arr, int* count, int index, size_t elem_size) {
    if (!arr || index < 0 || index >= *count) return -1;
    char* base = (char*)arr;
    if (index < *count - 1) {
        memmove(base + (size_t)index * elem_size,
                base + (size_t)(index + 1) * elem_size,
                (size_t)(*count - index - 1) * elem_size);
    }
    (*count)--;
    return 0;
}

static Joint make_joint(int id, double x, double y, int is_support, const char* releases) {
    Joint j = {0};
    j.id = id;
    j.x = x;
    j.y = y;
    j.is_support = is_support ? 1 : 0;
    if (releases) {
        snprintf(j.releases, sizeof(j.releases), "%s", releases);
    }
    return j;
}

static Member make_member(int id, int j1_id, int j2_id, double az, double iy, double ix) {
    return (Member){.id = id, .j1_id = j1_id, .j2_id = j2_id, .az = az, .iy = iy, .ix = ix};
}

static JointLoad make_joint_load(int jid, char axis, double magnitude) {
    return (JointLoad){.jid = jid, .axis = (char)toupper((unsigned char)axis), .magnitude = magnitude};
}

static MemberLoad make_member_load(int mid, char axis, char type, const double vals[4]) {
    MemberLoad ml = {0};
    ml.mid = mid;
    ml.axis = (char)toupper((unsigned char)axis);
    ml.type = (char)toupper((unsigned char)type);
    if (vals) {
        memcpy(ml.vals, vals, sizeof(ml.vals));
    }
    return ml;
}

/* StressModel setters & item mutators */
int stress_model_set_E(StressModel* model, double E) {
    if (!model) return -1;
    model->E = E;
    return 0;
}

int stress_model_set_G(StressModel* model, double G) {
    if (!model) return -1;
    model->G = G;
    return 0;
}

int stress_model_set_type(StressModel* model, const char* type) {
    if (!model || !type) return -1;
    snprintf(model->type, sizeof(model->type), "%s", type);
    return 0;
}

int stress_model_add_joint(StressModel* model, int id, double x, double y,
                           int is_support, const char* releases) {
    if (!model) return -1;
    return stress_model_insert_joint(model, model->num_joints, id, x, y, is_support, releases);
}

int stress_model_insert_joint(StressModel* model, int index, int id, double x, double y,
                              int is_support, const char* releases) {
    if (!model) return -1;
    Joint j = make_joint(id, x, y, is_support, releases);
    return insert_item((void**)&model->joints, &model->num_joints, &model->max_joints,
                       index, &j, sizeof(Joint));
}

int stress_model_update_joint(StressModel* model, int index, int id, double x, double y,
                              int is_support, const char* releases) {
    if (!model || index < 0 || index >= model->num_joints) return -1;
    model->joints[index] = make_joint(id, x, y, is_support, releases);
    return 0;
}

int stress_model_delete_joint(StressModel* model, int index) {
    if (!model) return -1;
    return delete_item(model->joints, &model->num_joints, index, sizeof(Joint));
}

int stress_model_add_member(StressModel* model, int id, int j1_id, int j2_id,
                            double az, double iy, double ix) {
    if (!model) return -1;
    return stress_model_insert_member(model, model->num_members, id, j1_id, j2_id, az, iy, ix);
}

int stress_model_insert_member(StressModel* model, int index, int id, int j1_id, int j2_id,
                               double az, double iy, double ix) {
    if (!model) return -1;
    Member m = make_member(id, j1_id, j2_id, az, iy, ix);
    return insert_item((void**)&model->members, &model->num_members, &model->max_members,
                       index, &m, sizeof(Member));
}

int stress_model_update_member(StressModel* model, int index, int id, int j1_id, int j2_id,
                               double az, double iy, double ix) {
    if (!model || index < 0 || index >= model->num_members) return -1;
    model->members[index] = make_member(id, j1_id, j2_id, az, iy, ix);
    return 0;
}

int stress_model_delete_member(StressModel* model, int index) {
    if (!model) return -1;
    return delete_item(model->members, &model->num_members, index, sizeof(Member));
}

int stress_model_add_joint_load(StressModel* model, int jid, char axis, double magnitude) {
    if (!model) return -1;
    return stress_model_insert_joint_load(model, model->num_joint_loads, jid, axis, magnitude);
}

int stress_model_insert_joint_load(StressModel* model, int index, int jid, char axis,
                                   double magnitude) {
    if (!model) return -1;
    JointLoad jl = make_joint_load(jid, axis, magnitude);
    return insert_item((void**)&model->joint_loads, &model->num_joint_loads,
                       &model->max_joint_loads, index, &jl, sizeof(JointLoad));
}

int stress_model_update_joint_load(StressModel* model, int index, int jid, char axis,
                                   double magnitude) {
    if (!model || index < 0 || index >= model->num_joint_loads) return -1;
    model->joint_loads[index] = make_joint_load(jid, axis, magnitude);
    return 0;
}

int stress_model_delete_joint_load(StressModel* model, int index) {
    if (!model) return -1;
    return delete_item(model->joint_loads, &model->num_joint_loads, index, sizeof(JointLoad));
}

int stress_model_add_member_load(StressModel* model, int mid, char axis, char type,
                                 const double vals[4]) {
    if (!model) return -1;
    return stress_model_insert_member_load(model, model->num_member_loads, mid, axis, type, vals);
}

int stress_model_insert_member_load(StressModel* model, int index, int mid, char axis,
                                    char type, const double vals[4]) {
    if (!model) return -1;
    MemberLoad ml = make_member_load(mid, axis, type, vals);
    return insert_item((void**)&model->member_loads, &model->num_member_loads,
                       &model->max_member_loads, index, &ml, sizeof(MemberLoad));
}

int stress_model_update_member_load(StressModel* model, int index, int mid, char axis,
                                    char type, const double vals[4]) {
    if (!model || index < 0 || index >= model->num_member_loads) return -1;
    model->member_loads[index] = make_member_load(mid, axis, type, vals);
    return 0;
}

int stress_model_delete_member_load(StressModel* model, int index) {
    if (!model) return -1;
    return delete_item(model->member_loads, &model->num_member_loads, index, sizeof(MemberLoad));
}
