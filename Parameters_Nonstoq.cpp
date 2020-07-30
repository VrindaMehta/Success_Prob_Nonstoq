#include "basic.h"
#include <iostream>
#include <string.h>

void params(string alpha , float g) {

    int i , j;

    for (i = 0; i < n; i++) {
            h_x[i] = 1;
        }

    for(i = 0; i < n; i++) {
         h1_x[i] = 0;
         h_y[i] = 0;
    }

    /*for(i = 0; i < n; i++) {   //To make the problem Hamiltonian nonstoquastic
      for(j = 0; j < n; j++) {
        if(flag[j + i * n] != 0) {
          J_y[j + i * n] = 0.5;
        }
      }
    } */

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
                if(flag[j + i * n] != 0) {
                    if (strcmp(alpha.c_str(),"A") == 0) {
                        J_x[j + i * n] = -g;
                    }
                    else if (strcmp(alpha.c_str(),"F") == 0) {
                        J_x[j + i * n] = g;
                    }
                    else if (strcmp(alpha.c_str(),"M") == 0) {
                        J_x[j + i * n] = float(rand() % 200 - 100) * g / 100;
                    }
                    else if (strcmp(alpha.c_str(),"O") == 0) {
                        J_x[j + i * n] = 0;
                    }
                }

                else {
                    J_x[j + i *n] = 0;
                }
        }
    }
}
