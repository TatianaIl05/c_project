#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include "clustering.h"

//finds the power two of the sum of squared differences between two points' coordinates
double euclidean_distance_squared(point a, point b, int dim)
{
    double res = 0;
    for (int i = 0; i < dim; i++)
    {
        res += (b.axis[i] - a.axis[i])*(b.axis[i] - a.axis[i]);
    }

    return res;
} 

//initialises points
void init_points(point * pnts, int N, int dim)
{
    for (int i = 0; i < N; i++)
    {
        pnts[i].axis = (double *)malloc(dim * sizeof(double));
    }
}

//frees the points
void free_points(point * pnts, int N)
{
    for (int i = 0; i < N; i++)
    {
        free(pnts[i].axis);
    }
}

/*
initialises centroids using the k-means++ algorithm:
randomly initialises the first centroid,
finds the square of the distance between each point and nearest centroid,
chooses the next centroid among these points so that the probability of its choice is proportional to the calculated squared distance,
finds all of the next centroids the same way
*/
void init_centroids(point * pnts, int N, point * centroids, int clust_num, int dim)
{
    int * indexes_for_centr = (int *)malloc(clust_num * sizeof(int));
    indexes_for_centr[0] = rand() % N;

    init_points(centroids, clust_num, dim); 
    for (int i = 0; i < dim; i++)
    {
        centroids[0].axis[i] = pnts[indexes_for_centr[0]].axis[i];
    }

    double * dst = (double *)malloc(N * sizeof(double));

    for (int i = 1; i < clust_num; i++)
    {
        double distance = 0;
        for (int j = 0; j < N; j++)
        {
            double dst_min = euclidean_distance_squared(centroids[0], pnts[j], dim);
            for (int k = 1; k < i; k++)
            {
                double new_dist = euclidean_distance_squared(pnts[j], centroids[k], dim);

                dst_min = (dst_min < new_dist) ? dst_min : new_dist;
            }

            distance += dst_min;
            dst[j] = dst_min;
        }

        int ind = 0;
        double rnd_num = distance * (rand() / (double) RAND_MAX);
        double break_point = 0;

        while (ind < N)
        {
            break_point += dst[ind++];
            if (break_point >= rnd_num) break;
        }
        --ind;
        indexes_for_centr[i] = ind;

        for (int j = 0; j < dim; j++)
        {
            centroids[i].axis[j] = pnts[ind].axis[j];
        }
    }

    free(indexes_for_centr);
    free(dst);
}

//takes two sorted arrays and combines them into one sorted array
void merge(point * pnts, int left, int middle, int right, int ind) 
{
    int first = middle - left + 1;
    int second = right - middle;

    point * leftArr = (point *)malloc(first * sizeof(point));
    point * rightArr = (point *)malloc(second * sizeof(point));

    for (int i = 0; i < first; i++)
    {
        leftArr[i] = pnts[left + i];
    }
    for (int j = 0; j < second; j++)
    {
        rightArr[j] = pnts[middle + 1 + j];
    }

    int i = 0, j = 0, k = left;

    while (i < first && j < second) 
    {
        if (leftArr[i].axis[ind] <= rightArr[j].axis[ind]) 
        {
            pnts[k] = leftArr[i++];
        }
        else 
        {
            pnts[k] = rightArr[j++];
        }
        k++;
    }

    while (i < first) 
    {
        pnts[k++] = leftArr[i++];
    }

    while (j < second) 
    {
        pnts[k++] = rightArr[j++];
    }

    free(leftArr);
    free(rightArr);
}

//divides an array into subarrays, sorts them and merges back together
void mergeSort(point * pnts, int left, int right, int ind) 
{
    if (left < right) 
    {
        int middle = (right - left) / 2 + left;

        mergeSort(pnts, left, middle, ind);
        mergeSort(pnts, middle + 1, right, ind);

        merge(pnts, left, middle, right, ind);
    }
}

//sorts the points based on ind and assignes the middle point to the root, recursively builds left and right subtrees
Node * make_kd_tree(point * pnts, int N, int * all_indexes, int ind, int dim)
{
    if (N <= 0) return NULL;

    mergeSort(pnts, 0, N - 1, ind % dim);

    int mid = N / 2;
    Node * point = (Node *)malloc(sizeof(Node));
    point->index = all_indexes[mid];
    point->coordinates.axis = (double *)malloc(dim * sizeof(double));
    
    for (int i = 0; i < dim; i++)
    {
        point->coordinates.axis[i] = pnts[mid].axis[i];
    }

    point->left = make_kd_tree(pnts, mid, all_indexes, (ind + 1) % dim, dim);
    point->right = make_kd_tree(pnts + 1 + mid, N - (mid + 1), all_indexes + 1 + mid, (ind + 1) % dim, dim);

    return point;
}

//while the current node is not NULL, recursively goes to the left and right subtree, frees the coordinates and the node
void kd_tree_free(Node * pnt)
{
    if (pnt == NULL) return;

    kd_tree_free(pnt->left);
    kd_tree_free(pnt->right);

    free(pnt->coordinates.axis);
    free(pnt);
}

/*
processes the point by comparing it with the current node of the kd-tree, if the distance between them is less than the best distance so far, 
changes the best distance to the current one, considering the necessary dimention compares the point in the node and the 
processing point, if the processing point is smaller, recursively continues the same in the left subtree, otherwise continues 
in the right one
*/
void nearest_neighbour(Node * pnt, point processing, double * dst, int ind, int dim, int * best_ind)
{
    if (pnt == NULL) return; 

    double dist = euclidean_distance_squared(pnt->coordinates, processing, dim);

    if (dist < *dst)
    {
        *dst = dist;
        *best_ind = pnt->index;
    } 

    if (processing.axis[ind % dim] < pnt->coordinates.axis[ind % dim])
    {
        nearest_neighbour(pnt->left, processing, dst, ind + 1, dim, best_ind);

        if (abs(processing.axis[ind % dim] - pnt->coordinates.axis[ind % dim]) < *dst)
        {
            nearest_neighbour(pnt->right, processing, dst, ind + 1, dim, best_ind);
        }
    }
    else
    {
        nearest_neighbour(pnt->right, processing, dst, ind + 1, dim, best_ind);
        
        if (abs(processing.axis[ind % dim] - pnt->coordinates.axis[ind % dim]) < *dst)
        {
            nearest_neighbour(pnt->left, processing, dst, ind + 1, dim, best_ind);
        }
    }
}

/*
intialises centroids, makes a kd-tree, for each point finds the nearest centroid, if the centroids are the same as the ones from 
the previous iteration, it stops, otherwise it updates the centroids by calculating the centre of mass or choosing a new random centroid 
if there are no points related to the previous one
*/
void k_means(point * pnts, int N, point * centroids, int clust_num, int dim, int iterat_num, int * indexes, int * indec)
{
    init_centroids(pnts, N, centroids, clust_num, dim);
    
    for (int i = 0; i < iterat_num; i++)
    {
        int * indexes = (int *)malloc(clust_num * sizeof(int));
        for (int j = 0; j < clust_num; j++)
        {
            indexes[j] = j;
        }

        int stop = 1;

        Node * root = make_kd_tree(centroids, clust_num, indexes, 0, dim);

        for (int j = 0; j < N; j++)
        {
            point processing = pnts[j];
            double dist = DBL_MAX;
            int best_ind = -1;
    
            nearest_neighbour(root, processing, &dist, 0, dim, &best_ind);

            if (indec[j] != best_ind) stop = 0;
            indec[j] = best_ind;
        }

        free(indexes);
        kd_tree_free(root);
        
        if (stop) break;

        for (int j = 0; j < clust_num; j++)
        {
            double * sums = (double *)calloc(dim, sizeof(double));
            int change = 0;
            for (int k = 0; k < N; k++)
            {
                if (indec[k] == j)
                {
                    for (int l = 0; l < dim; l++)
                    {
                        sums[l] += pnts[k].axis[l];
                        
                    }
                    change++;
                }
            }
            
            if (change) 
            {
                for (int k = 0; k < dim; k++)
                {
                    centroids[j].axis[k] = sums[k] / (double)change;
                }
            }
            else
            {
                int random = rand() % N;
                for (int k = 0; k < dim; k++)
                {
                    centroids[j].axis[k] = pnts[random].axis[k];
                }
            }

            free(sums);
        }
        
    }
}

/*
processes the point by comparing it with the current node of the kd-tree, if the distance between them is less than epsilon, 
includes the point in the node to the neighbours, considering the necessary dimention compares the point in the node and the 
processing point, if the processing point is smaller, recursively continues the same in the left subtree, otherwise continues 
in the right one, backtracks
*/
void all_neighbours(Node * pnt, point processing, int * neigh_count, int * indexes, double epsilon, int ind, int dim)
{
    if (pnt == NULL) return;

    double distance = euclidean_distance_squared(pnt->coordinates, processing, dim);
    
    if (distance <= epsilon)
    {
        indexes[*neigh_count] = pnt->index;
        (*neigh_count)++;
    }
    
    if (processing.axis[ind] <= pnt->coordinates.axis[ind])
    {
        all_neighbours(pnt->left, processing, neigh_count, indexes, epsilon, (ind + 1) % dim, dim);

        if (abs(processing.axis[ind % dim] - pnt->coordinates.axis[ind % dim]) <= epsilon)
        {
            all_neighbours(pnt->right, processing, neigh_count, indexes, epsilon, (ind + 1) % dim, dim);
        }
    }
    else
    {
        all_neighbours(pnt->right, processing, neigh_count, indexes, epsilon, (ind + 1) % dim, dim);
        
        if (abs(processing.axis[ind % dim] - pnt->coordinates.axis[ind % dim]) <= epsilon)
        {
            all_neighbours(pnt->left, processing, neigh_count, indexes, epsilon, (ind + 1) % dim, dim);
        }
    }

}

/*
assigns the current point to the current cluster, goes through all of current point's neighbours and finds their neighbours, if the number of their
neighbours makes them a main point, continues making the same cluster recursively, assigns those neighbours to the current cluster
*/
void making_cluster(Node * root, point * pnts, int N, int * neigh_count, int * indexes, double epsilon, int min_density, int * visit, int dim, int * indec, int clust_num, int current_ind)
{
    indec[current_ind] = clust_num;
    for (int i = 0; i < *neigh_count; i++)
    {
        if (!visit[indexes[i]])
        {
            visit[indexes[i]] = 1;

            int * new_indexes = (int *)malloc(N * sizeof(int));
            int new_neigh_count = 0;

            all_neighbours(root, pnts[indexes[i]], &new_neigh_count, new_indexes, epsilon, 0, dim);

            if (new_neigh_count >= min_density)
            {
                making_cluster(root, pnts, N, &new_neigh_count, new_indexes, epsilon, min_density, visit, dim, indec, clust_num, indexes[i]);
            }

            free(new_indexes);
        }

        if (indec[indexes[i]] == 0) indec[indexes[i]] = clust_num;
    }
}

/*
visits one by one aech of the unvisited points, finds all its neighbours and if the point is main, makes a cluster with it,
otherwise marks it as an outlier
*/
int dbscan(point * pnts, int N, double epsilon, int min_density, int dim, int * indec)
{
    int * indexes_in_mass = (int *)malloc(N * sizeof(int));
    for (int i = 0; i < N; i++)
    {
        indexes_in_mass[i] = i;
    }

    int clust_num = 1;
    int * visit = (int *)calloc(N, sizeof(int));
    Node * root = make_kd_tree(pnts, N, indexes_in_mass, 0, dim);

    for (int i = 0; i < N; i++)
    {
        if (!visit[i])
        {
            visit[i] = 1;

            int * indexes_in_init = (int *)malloc(N * sizeof(int));
            int neigh_count = 0;
            all_neighbours(root, pnts[i], &neigh_count, indexes_in_init, epsilon, 0, dim);

            if (neigh_count < min_density)
            {
                indec[i] = -1;
            }
            else
            {
                indec[i] = clust_num;
                making_cluster(root, pnts, N, &neigh_count, indexes_in_init, epsilon, min_density, visit, dim, indec, clust_num, i);
                clust_num++;
            }

            free(indexes_in_init);
        }
    }

    kd_tree_free(root);
    free(indexes_in_mass);
    free(visit);

    return clust_num;
}

/*
processes each point of an array pnts by considering its cluster and adding the appropriate coordinate 
to its sum, divides the sum by the number of points which belong to a certain cluster 
*/
void clusters_for_dbscan(point ** centres, point * pnts, int clust_count, int N, int dim, int * indec)
{
    *centres = (point *)malloc(clust_count * sizeof(point));
    init_points(*centres, clust_count, dim);
    for (int i = 0; i < clust_count; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            (*centres)[i].axis[j] = 0;
        }
    }

    int * counts = (int *)calloc(clust_count, sizeof(int));
    for (int i = 0; i < N; i++)
    {
        if (indec[i] != -1)
        {
            for (int j = 0; j < dim; j++)
            {
                (*centres)[indec[i]].axis[j] += pnts[i].axis[j];
            }
            counts[indec[i]]++;
        }
    }

    for (int i = 1; i < clust_count; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            (*centres)[i].axis[j] /= counts[i];
        }
    }  

    free(counts);
}

/*
reads the first row (it should be less than 1024 characters long), counts the number of axes by considering 
the number of spaces between words in the first row (strictly one space between each word) and assuming this 
dim stores the data into an array
*/
point * data_set_read(int N, int * dim)
{
    FILE * input = fopen("input.txt", "r");

    char space;
    char ch[1024];
    *dim = 0;
    if (fgets(ch, sizeof(ch), input) != NULL)
    {
        for (int i = 0; ch[i] != '\0'; i++)
        {
            if (ch[i] == ' ') (*dim)++;
        }
        (*dim)++;
    } 

    point * data = (point *)malloc(N * sizeof(point));
    init_points(data, N, *dim);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < *dim; j++ )
        {
            fscanf(input, "%lf", &data[i].axis[j]);
            fscanf(input, "%c", &space);
        }
    }
    
    fclose(input);

    return data;
}

//prints all of the centroids
void printing_centroids(point * centroids, int clust_num, int dim, int flag)
{
    printf("Centroids:\n");
    
    for (int i = flag; i < clust_num; i++)
    {
        printf("Centroid %d - (", i-flag);
        for (int j = 0; j < dim; j++)
        {
            if (j == dim - 1) printf("%.3f)\n", centroids[i].axis[j]);
            else printf("%.3f, ", centroids[i].axis[j]);
        }
    }
}

//prints all of the points along with the number of their cluster or the flag "Outlier"
void printing_points(point * points, int N, int dim, int * indec, int flag)
{
    for (int i = 0; i < N; i++)
    {
        printf("Point (");

        for (int j = 0; j < dim; j++)
        {
            if (j == dim - 1) printf("%.3f) -> ", points[i].axis[j]);
            else printf("%.3f, ", points[i].axis[j]);
        }

        if (indec[i] == -1) printf("Outlier\n");
        else printf("Cluster %d\n", indec[i] - flag);
    }
}