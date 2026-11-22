#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Joint Joint;
typedef struct Member Member;
typedef struct JointLoad JointLoad;
typedef struct MemberLoad MemberLoad;
typedef struct StressModel StressModel;

StressModel* stress_create(void);
int stress_load_dat(StressModel* model, const char* filename);
int stress_solve(StressModel* model);
void stress_free(StressModel* model);

/* Joint getters */
int stress_joint_get_id(const Joint* joint);
double stress_joint_get_x(const Joint* joint);
double stress_joint_get_y(const Joint* joint);
int stress_joint_get_is_support(const Joint* joint);
const char* stress_joint_get_releases(const Joint* joint);

/* Member getters */
int stress_member_get_id(const Member* member);
int stress_member_get_j1_id(const Member* member);
int stress_member_get_j2_id(const Member* member);
double stress_member_get_az(const Member* member);
double stress_member_get_iy(const Member* member);
double stress_member_get_ix(const Member* member);

/* JointLoad getters */
int stress_joint_load_get_jid(const JointLoad* load);
char stress_joint_load_get_axis(const JointLoad* load);
double stress_joint_load_get_magnitude(const JointLoad* load);

/* MemberLoad getters */
int stress_member_load_get_mid(const MemberLoad* load);
char stress_member_load_get_axis(const MemberLoad* load);
char stress_member_load_get_type(const MemberLoad* load);
const double* stress_member_load_get_vals(const MemberLoad* load);

/* StressModel getters */
const Joint* stress_model_get_joint(const StressModel* model, int index);
int stress_model_get_num_joints(const StressModel* model);
int stress_model_get_max_joints(const StressModel* model);

const Member* stress_model_get_member(const StressModel* model, int index);
int stress_model_get_num_members(const StressModel* model);
int stress_model_get_max_members(const StressModel* model);

const JointLoad* stress_model_get_joint_load(const StressModel* model, int index);
int stress_model_get_num_joint_loads(const StressModel* model);
int stress_model_get_max_joint_loads(const StressModel* model);

const MemberLoad* stress_model_get_member_load(const StressModel* model, int index);
int stress_model_get_num_member_loads(const StressModel* model);
int stress_model_get_max_member_loads(const StressModel* model);

double stress_model_get_E(const StressModel* model);
double stress_model_get_G(const StressModel* model);
const char* stress_model_get_type(const StressModel* model);

const double* stress_model_get_displacements(const StressModel* model);
const double* stress_model_get_reactions(const StressModel* model);

#ifdef __cplusplus
}
#endif
