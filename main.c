#include "stress.h"
#include <stdio.h>

static void stress_print_results(const StressModel* model) {
    int num_joints = stress_model_get_num_joints(model);
    const double* displacements = stress_model_get_displacements(model);
    const double* reactions = stress_model_get_reactions(model);

    printf("\nJOINT DISPLACEMENTS (Z, RotX, RotY):\n");
    printf("%-8s %12s %12s %12s\n", "Joint", "Z-Disp", "Rot-X", "Rot-Y");
    for (int i = 0; i < num_joints; i++) {
        const Joint* j = stress_model_get_joint(model, i);
        printf("%-8d %12.6e %12.6e %12.6e\n", stress_joint_get_id(j),
               displacements[i*3], displacements[i*3+1], displacements[i*3+2]);
    }
    printf("\nSUPPORT REACTIONS:\n");
    printf("%-8s %12s %12s %12s\n", "Joint", "F-Z", "M-X", "M-Y");
    for (int i = 0; i < num_joints; i++) {
        const Joint* j = stress_model_get_joint(model, i);
        if (stress_joint_get_is_support(j)) {
            printf("%-8d %12.2f %12.2f %12.2f\n", stress_joint_get_id(j),
                   reactions[i*3], reactions[i*3+1], reactions[i*3+2]);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s <filename.dat>\n", argv[0]);
        return 1;
    }

    StressModel* model = stress_create();
    if (!model) {
        printf("Error allocating model.\n");
        return 1;
    }
    if (stress_load_dat(model, argv[1]) != 0) {
        printf("Error loading file: %s\n", argv[1]);
        stress_free(model);
        return 1;
    }

    printf("Solving %s model with %d joints and %d members.\n",
           stress_model_get_type(model),
           stress_model_get_num_joints(model),
           stress_model_get_num_members(model));

    if (stress_solve(model) == 0) {
        stress_print_results(model);
    } else {
        printf("Error solving model.\n");
    }

    stress_free(model);
    return 0;
}
