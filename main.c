#include "stress.h"
#include <stdio.h>

void stress_print_results(StressModel* model) {
    printf("\nJOINT DISPLACEMENTS (Z, RotX, RotY):\n");
    printf("%-8s %12s %12s %12s\n", "Joint", "Z-Disp", "Rot-X", "Rot-Y");
    for (int i = 0; i < model->num_joints; i++) {
        printf("%-8d %12.6e %12.6e %12.6e\n", model->joints[i].id, 
               model->displacements[i*3], model->displacements[i*3+1], model->displacements[i*3+2]);
    }
    printf("\nSUPPORT REACTIONS:\n");
    printf("%-8s %12s %12s %12s\n", "Joint", "F-Z", "M-X", "M-Y");
    for (int i = 0; i < model->num_joints; i++) {
        if (model->joints[i].is_support) {
            printf("%-8d %12.2f %12.2f %12.2f\n", model->joints[i].id, 
                   model->reactions[i*3], model->reactions[i*3+1], model->reactions[i*3+2]);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Usage: %s <filename.dat>\n", argv[0]);
        return 1;
    }

    StressModel* model = stress_create();
    if (stress_load_dat(model, argv[1]) != 0) {
        printf("Error loading file: %s\n", argv[1]);
        stress_free(model);
        return 1;
    }

    printf("Solving %s model with %d joints and %d members.\n", 
           model->type, model->num_joints, model->num_members);

    if (stress_solve(model) == 0) {
        stress_print_results(model);
    } else {
        printf("Error solving model.\n");
    }

    stress_free(model);
    return 0;
}
