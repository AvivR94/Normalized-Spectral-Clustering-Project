#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

int main() {
    double w_matrix[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
    double vectors_matrix[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
    //calculating the values in each cell
    int sum = 0;
    for (int i = 0; i < 3; i++)
        w_matrix[i][i] = 0; //the diagonal line in matrix's values are 0
    for (int i = 0; i < 3; i++){
        for (int j= i + 1; j < 3; j++){
            sum = 0;
            for (int s = 0; s < 3; s++)
                sum += pow((vectors_matrix[i][s] - vectors_matrix[j][s]),2);
            w_matrix[i][j] = exp(pow(sum,(0.5))/-2);
            w_matrix[j][i] = exp(pow(sum,(0.5))/-2);
        }
    }
    printf("fgd");
    return 0;
}