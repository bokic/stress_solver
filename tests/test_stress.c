#include "stress.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(cond)                                                             \
    do {                                                                        \
        if (!(cond)) {                                                          \
            fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__,  \
                    #cond);                                                     \
            exit(1);                                                            \
        }                                                                       \
    } while (0)

static int approx(double a, double b, double tol) {
    return fabs(a - b) <= tol;
}

static StressModel* load_from_string(const char* content) {
    const char* path = "/tmp/stress_test_tmp.dat";
    FILE* f = fopen(path, "w");
    CHECK(f != NULL);
    fputs(content, f);
    fclose(f);

    StressModel* m = stress_create();
    CHECK(m != NULL);
    CHECK(stress_load_dat(m, path) == 0);
    remove(path);
    return m;
}

static void test_cantilever_point_load(void) {
    StressModel* m = load_from_string(
        "TYPE PLANE GRID\n"
        "CONSTANTS E 2.0e7\n"
        "JOINT COOR\n"
        "1 0.0 0.0 SUPPORT\n"
        "2 4.0 0.0\n"
        "MEMB INCI\n"
        "1 1 2\n"
        "MEMBER PROPERTIES\n"
        "1 AZ 1.0 IY 1.0e-3 IX 1.0e-3\n"
        "JOINT LOADS\n"
        "2 FORCE Z -10.0\n"
        "SOLVE\n");

    CHECK(stress_solve(m) == 0);

    double E = stress_model_get_E(m);
    const double* disp = stress_model_get_displacements(m);
    const double* react = stress_model_get_reactions(m);

    // w = P*L^3 / (3*E*I)
    double expected_w = -10.0 * 64.0 / (3.0 * E * 1.0e-3);
    CHECK(approx(disp[3], expected_w, 1e-9));
    // Vertical reaction at support 1 = +10.0
    CHECK(approx(react[0], 10.0, 1e-6));

    // Re-solve should succeed cleanly without leaks
    CHECK(stress_solve(m) == 0);
    disp = stress_model_get_displacements(m);
    CHECK(approx(disp[3], expected_w, 1e-9));

    stress_free(m);
}

static void test_simply_supported_uniform_load(void) {
    StressModel* m = load_from_string(
        "TYPE PLANE GRID\n"
        "CONSTANTS E 1.0e7\n"
        "JOINT COOR\n"
        "1 0.0 0.0 SUPPORT\n"
        "2 3.0 0.0\n"
        "3 6.0 0.0 SUPPORT\n"
        "MEMB INCI\n"
        "1 1 2\n"
        "2 2 3\n"
        "MEMBER PROPERTIES\n"
        "1 THRU 2 AZ 1.0 IY 1.0e-3 IX 1.0e-3\n"
        "JOINT RELEASES\n"
        "1 MOMENT Y\n"
        "3 MOMENT Y\n"
        "MEMBER LOADS\n"
        "1 THRU 2 FORCE Z UNIFORM   -2.0\n"
        "SOLVE\n");

    CHECK(strcmp(stress_model_get_type(m), "PLANE GRID") == 0);
    CHECK(strcmp(stress_joint_get_releases(stress_model_get_joint(m, 0)), "MY") == 0);
    CHECK(stress_solve(m) == 0);

    const double* disp = stress_model_get_displacements(m);
    const double* react = stress_model_get_reactions(m);

    // Total load = -2.0 * 6.0 = -12.0 -> reactions at 1 and 3 are +6.0 each
    CHECK(approx(react[0], 6.0, 1e-5));
    CHECK(approx(react[6], 6.0, 1e-5));

    // Midspan deflection: 5 * w * L^4 / (384 * E * I)
    double L = 6.0;
    double expected_mid = 5.0 * (-2.0) * pow(L, 4) / (384.0 * stress_model_get_E(m) * 1.0e-3);
    CHECK(approx(disp[3], expected_mid, 1e-8));

    stress_free(m);
}

static void test_getters_and_dynamic_growth(void) {
    const char* path = "/tmp/stress_test_large.dat";
    FILE* f = fopen(path, "w");
    CHECK(f != NULL);
    fprintf(f, "TYPE PLANE GRID\nCONSTANTS E 2.0e7 POISSON 0.25\nJOINT COOR\n");
    for (int i = 1; i <= 120; i++) {
        fprintf(f, "%d %.2f 0.0%s\n", i, (double)(i - 1), (i == 1) ? " SUPPORT" : "");
    }
    fprintf(f, "MEMB INCI\n");
    for (int i = 1; i <= 119; i++) {
        fprintf(f, "%d %d %d\n", i, i, i + 1);
    }
    fprintf(f, "MEMBER PROPERTIES\n1 THRU 119 AZ 1.0 IY 1.0e-3 IX 2.0e-3\n");
    fprintf(f, "JOINT LOADS\n120 FORCE Z -1.0\n");
    fprintf(f, "MEMBER LOADS\n1 FORCE Z CONCENTRATED -2.0 0.5\nSOLVE\n");
    fclose(f);

    StressModel* m = stress_create();
    CHECK(stress_load_dat(m, path) == 0);
    remove(path);

    CHECK(stress_model_get_num_joints(m) == 120);
    CHECK(stress_model_get_max_joints(m) >= 120);
    CHECK(stress_model_get_num_members(m) == 119);
    CHECK(stress_model_get_max_members(m) >= 119);
    CHECK(stress_model_get_num_joint_loads(m) == 1);
    CHECK(stress_model_get_max_joint_loads(m) >= 1);
    CHECK(stress_model_get_num_member_loads(m) == 1);
    CHECK(stress_model_get_max_member_loads(m) >= 1);
    CHECK(approx(stress_model_get_E(m), 2.0e7, 1e-6));
    CHECK(approx(stress_model_get_G(m), 2.0e7 / 2.5, 1e-3));

    const Joint* j1 = stress_model_get_joint(m, 0);
    CHECK(stress_joint_get_id(j1) == 1);
    CHECK(approx(stress_joint_get_x(j1), 0.0, 1e-9));
    CHECK(approx(stress_joint_get_y(j1), 0.0, 1e-9));
    CHECK(stress_joint_get_is_support(j1) == 1);
    CHECK(strcmp(stress_joint_get_releases(j1), "") == 0);

    const Member* m1 = stress_model_get_member(m, 0);
    CHECK(stress_member_get_id(m1) == 1);
    CHECK(stress_member_get_j1_id(m1) == 1);
    CHECK(stress_member_get_j2_id(m1) == 2);
    CHECK(approx(stress_member_get_az(m1), 1.0, 1e-9));
    CHECK(approx(stress_member_get_iy(m1), 1.0e-3, 1e-9));
    CHECK(approx(stress_member_get_ix(m1), 2.0e-3, 1e-9));

    const JointLoad* jl = stress_model_get_joint_load(m, 0);
    CHECK(stress_joint_load_get_jid(jl) == 120);
    CHECK(stress_joint_load_get_axis(jl) == 'Z');
    CHECK(approx(stress_joint_load_get_magnitude(jl), -1.0, 1e-9));

    const MemberLoad* ml = stress_model_get_member_load(m, 0);
    CHECK(stress_member_load_get_mid(ml) == 1);
    CHECK(stress_member_load_get_axis(ml) == 'Z');
    CHECK(stress_member_load_get_type(ml) == 'C');
    const double* vals = stress_member_load_get_vals(ml);
    CHECK(approx(vals[0], -2.0, 1e-9));
    CHECK(approx(vals[1], 0.5, 1e-9));

    CHECK(stress_solve(m) == 0);
    CHECK(approx(stress_model_get_reactions(m)[0], 3.0, 1e-3));
    stress_free(m);
}

static void test_error_cases(void) {
    // Singular matrix (no supports)
    StressModel* sing = load_from_string(
        "TYPE PLANE GRID\nCONSTANTS E 2.0e7\n"
        "JOINT COOR\n1 0.0 0.0\n2 1.0 0.0\n"
        "MEMB INCI\n1 1 2\n"
        "MEMBER PROPERTIES\n1 AZ 1.0 IY 1.0 IX 1.0\nSOLVE\n");
    CHECK(stress_solve(sing) == -1);
    stress_free(sing);

    // Zero-length member
    StressModel* zl = load_from_string(
        "TYPE PLANE GRID\nCONSTANTS E 2.0e7\n"
        "JOINT COOR\n1 0.0 0.0 SUPPORT\n2 0.0 0.0\n"
        "MEMB INCI\n1 1 2\n"
        "MEMBER PROPERTIES\n1 AZ 1.0 IY 1.0 IX 1.0\nSOLVE\n");
    CHECK(stress_solve(zl) == -1);
    stress_free(zl);

    // Missing joint ID
    StressModel* mj = load_from_string(
        "TYPE PLANE GRID\nCONSTANTS E 2.0e7\n"
        "JOINT COOR\n1 0.0 0.0 SUPPORT\n"
        "MEMB INCI\n1 1 99\n"
        "MEMBER PROPERTIES\n1 AZ 1.0 IY 1.0 IX 1.0\nSOLVE\n");
    CHECK(stress_solve(mj) == -1);
    stress_free(mj);
}

int main(void) {
    test_cantilever_point_load();
    test_simply_supported_uniform_load();
    test_getters_and_dynamic_growth();
    test_error_cases();
    printf("All stress_solver tests passed.\n");
    return 0;
}
