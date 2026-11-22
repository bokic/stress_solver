#include "stress.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int approx(double a, double b, double tol) {
    return fabs(a - b) <= tol;
}

static void test_cantilever_point_load(void) {
    StressModel* m = stress_create();
    assert(m != NULL);
    m->E = 2.0e7;
    m->G = m->E / 2.4;

    m->joints[0] = (Joint){.id = 1, .x = 0.0, .y = 0.0, .is_support = 1};
    m->joints[1] = (Joint){.id = 2, .x = 4.0, .y = 0.0, .is_support = 0};
    m->num_joints = 2;

    m->members[0] = (Member){.id = 1, .j1_id = 1, .j2_id = 2, .az = 1.0, .iy = 1.0e-3, .ix = 1.0e-3};
    m->num_members = 1;

    m->joint_loads[0] = (JointLoad){.jid = 2, .axis = 'Z', .magnitude = -10.0};
    m->num_joint_loads = 1;

    assert(stress_solve(m) == 0);

    // w = P*L^3 / (3*E*I) = -10 * 64 / (3 * 2e7 * 1e-3) = -640 / 60000 = -0.0106666667
    double expected_w = -10.0 * 64.0 / (3.0 * m->E * 1.0e-3);
    assert(approx(m->displacements[3], expected_w, 1e-9));
    // Vertical reaction at support 1 = +10.0
    assert(approx(m->reactions[0], 10.0, 1e-6));

    // Re-solve should succeed cleanly without leaks
    assert(stress_solve(m) == 0);
    assert(approx(m->displacements[3], expected_w, 1e-9));

    stress_free(m);
}

static void test_simply_supported_uniform_load(void) {
    const char* path = "/tmp/stress_test_ss.dat";
    FILE* f = fopen(path, "w");
    assert(f != NULL);
    fprintf(f,
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
    fclose(f);

    StressModel* m = stress_create();
    assert(stress_load_dat(m, path) == 0);
    assert(strcmp(m->type, "PLANE GRID") == 0);
    assert(stress_solve(m) == 0);

    // Total load = -2.0 * 6.0 = -12.0 -> reactions at 1 and 3 are +6.0 each
    assert(approx(m->reactions[0], 6.0, 1e-5));
    assert(approx(m->reactions[6], 6.0, 1e-5));

    // Midspan deflection: 5 * w * L^4 / (384 * E * I)
    double L = 6.0;
    double expected_mid = 5.0 * (-2.0) * pow(L, 4) / (384.0 * m->E * 1.0e-3);
    assert(approx(m->displacements[3], expected_mid, 1e-8));

    stress_free(m);
    remove(path);
}

static void test_dynamic_growth_and_errors(void) {
    // > 100 joints and > 200 members
    const char* path = "/tmp/stress_test_large.dat";
    FILE* f = fopen(path, "w");
    assert(f != NULL);
    fprintf(f, "TYPE PLANE GRID\nCONSTANTS E 2.0e7\nJOINT COOR\n");
    for (int i = 1; i <= 120; i++) {
        fprintf(f, "%d %.2f 0.0%s\n", i, (double)(i - 1), (i == 1) ? " SUPPORT" : "");
    }
    fprintf(f, "MEMB INCI\n");
    for (int i = 1; i <= 119; i++) {
        fprintf(f, "%d %d %d\n", i, i, i + 1);
    }
    fprintf(f, "MEMBER PROPERTIES\n1 THRU 119 AZ 1.0 IY 1.0e-3 IX 1.0e-3\n");
    fprintf(f, "JOINT LOADS\n120 FORCE Z -1.0\nSOLVE\n");
    fclose(f);

    StressModel* m = stress_create();
    assert(stress_load_dat(m, path) == 0);
    assert(m->num_joints == 120);
    assert(stress_solve(m) == 0);
    assert(approx(m->reactions[0], 1.0, 1e-4));
    stress_free(m);
    remove(path);

    // Singular matrix (no supports)
    StressModel* sing = stress_create();
    sing->E = 2.0e7;
    sing->G = 1.0e7;
    sing->joints[0] = (Joint){.id = 1, .x = 0.0, .y = 0.0, .is_support = 0};
    sing->joints[1] = (Joint){.id = 2, .x = 1.0, .y = 0.0, .is_support = 0};
    sing->num_joints = 2;
    sing->members[0] = (Member){.id = 1, .j1_id = 1, .j2_id = 2, .az = 1.0, .iy = 1.0, .ix = 1.0};
    sing->num_members = 1;
    assert(stress_solve(sing) == -1);
    stress_free(sing);

    // Zero-length member
    StressModel* zl = stress_create();
    zl->E = 2.0e7;
    zl->G = 1.0e7;
    zl->joints[0] = (Joint){.id = 1, .x = 0.0, .y = 0.0, .is_support = 1};
    zl->joints[1] = (Joint){.id = 2, .x = 0.0, .y = 0.0, .is_support = 0};
    zl->num_joints = 2;
    zl->members[0] = (Member){.id = 1, .j1_id = 1, .j2_id = 2, .az = 1.0, .iy = 1.0, .ix = 1.0};
    zl->num_members = 1;
    assert(stress_solve(zl) == -1);
    stress_free(zl);

    // Missing joint ID
    StressModel* mj = stress_create();
    mj->E = 2.0e7;
    mj->G = 1.0e7;
    mj->joints[0] = (Joint){.id = 1, .x = 0.0, .y = 0.0, .is_support = 1};
    mj->num_joints = 1;
    mj->members[0] = (Member){.id = 1, .j1_id = 1, .j2_id = 99, .az = 1.0, .iy = 1.0, .ix = 1.0};
    mj->num_members = 1;
    assert(stress_solve(mj) == -1);
    stress_free(mj);
}

int main(void) {
    test_cantilever_point_load();
    test_simply_supported_uniform_load();
    test_dynamic_growth_and_errors();
    printf("All stress_solver tests passed.\n");
    return 0;
}
