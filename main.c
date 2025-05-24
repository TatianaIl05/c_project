#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "clustering.h"

int main()
{
    srand(time(NULL)); 

    int dim, n = 100, clust_num = 8;

    point * data = data_set_read(n, &dim);

    //k_means
    point * centroids = (point *)malloc(clust_num * sizeof(point));
    int * indexes;
    int *indec1 = (int *)calloc(n, sizeof(int));
    
    k_means(data, n, centroids, clust_num, dim, 100, indexes, indec1);

    printing_centroids(centroids, clust_num, dim, 0);
    printing_points(data, n, dim, indec1, 0);

    free_points(centroids, clust_num);
    free(centroids);
    free(indec1);
    
    printf("\n");
    
    //dbscan
    int min_density = 5;
    double epsilon = (10 * 1 * 5) * (10 * 1 * 5);
    int *indec2 = (int *)calloc(n, sizeof(int));
    point * centroids2;

    int clust_count_pl1 = dbscan(data, n, epsilon, min_density, dim, indec2);
    clusters_for_dbscan(&centroids2, data, clust_count_pl1, n, dim, indec2);

    printing_centroids(centroids2, clust_count_pl1, dim, 1);
    printing_points(data, n, dim, indec2, 1);

    free(indec2);
    free_points(centroids2, clust_count_pl1);
    free(centroids2);


   free_points(data, n);  
   free(data);   
}