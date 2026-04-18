#ifndef STRESS_LIB_H
#define STRESS_LIB_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int id;
    double x, y;
    int is_support;
    char releases[8]; // e.g. "MX", "MY"
} Joint;

typedef struct {
    int id;
    int j1_id, j2_id;
    double az, iy, ix;
} Member;

typedef struct {
    int jid;
    char axis; // 'Z', 'X', 'Y'
    double magnitude;
} JointLoad;

typedef struct {
    int mid;
    char axis;
    char type; // 'U' for Uniform, 'C' for Conc
    double vals[4];
} MemberLoad;

typedef struct {
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
} StressModel;

StressModel* stress_create();
int stress_load_dat(StressModel* model, const char* filename);
int stress_solve(StressModel* model);
void stress_free(StressModel* model);

#ifdef __cplusplus
}
#endif

#endif // STRESS_LIB_H
