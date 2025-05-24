// structure for a point
typedef struct {
    double * axis; //number of axes in space (number of considered arguments of a point)
}point;

//structure for a node in a kd-tree
typedef struct Node{
    int index; //index in the tree
    point coordinates; //point in the node of a tree
    struct Node * left, * right; //pointers to left and right subtree
} Node;

/*
overview: counts the power two of the euclidean distance 
          between points according to the formula sqrt((x2 - x1)^2 + (y2 - y1)^2 + ...)
arguments: a - first point, b - second point, dim - number of dimentions (axes)
return: the sum of squared axis coordinates
*/
double euclidean_distance_squared(point a, point b, int dim);

/*
overview: initialises points 
arguments: pnts - pointer to an array of points, N - size of the array, dim - number of axes
*/
void init_points(point * pnts, int N, int dim);

/*
overview: frees the memory allocated in the init_points function
arguments: pnts - pointer to an array of points, N - size of the array
*/
void free_points(point * pnts, int N);

/*
overview: initialises centroids for the beginning of the k-means algorithm
arguments: pnts - pointer to an array of points that need to be divided into clusters, 
           N - size of array pnts, centroids - pointer to an array of centroids,
           clust_num - number of clusters, dim - number of axes
*/
void init_centroids(point * pnts, int N, point * centroids, int clust_num, int dim);

/*
overview: merges two subarrays
arguments: pnts - pointer to an array of points, left - starting index of the left subarray,
           middle - ending index of the left and starting index of the right subarrays,
           right - ending index of the right subarray, 
           ind - index of the axis used as a sorting criteria
*/
void merge(point * pnts, int left, int middle, int right, int ind);

/*
overview: sorts an array of points using the merge sort algorithm to build a kd-tree
arguments: pnts - pointer to an array of points, left - starting index of the subarray being sorted,
           right - ending index of the subarray being sorted, 
           ind - index of the axis used as a sorting criteria
*/
void mergeSort(point * pnts, int left, int right, int ind);

/*
overview: creates a kd-tree
arguments: pnts - pointer to an array of points, N - size of array pnts, 
           all_indexes - pointer to an array of indexes of points, ind - index of the axis used as a sorting criteria 
           to create new tree level, dim - number of axes
return: root of the created tree
*/
Node * make_kd_tree(point * pnts, int N, int * all_indexes, int ind, int dim);

/*
overview: frees all nodes of a kd-tree
arguments: pnt - pointer to a root of a kd-tree
*/
void kd_tree_free(Node * pnt);

/*
overview: finds the nearest centroid to the point to include it in the cluster related to this centroid
arguments: pnt - pointer to the root of kd-tree of centroids, processing - point for which the nearest centriod is being found,
           dst - pointer to a double where the distance to the nearest centroid is being stored,
           ind - index of the dimention used for comparison, dim - number of axes, best_ind - pointer to an int where the index of 
           the nearest centroid is being stored
*/
void nearest_neighbour(Node * pnt, point processing, double * dst, int ind, int dim, int * best_ind);

/*
overview: performs a k_means algorithm dividing pnts into clust_num clusters
arguments: pnts - pointer to an array of points, N - size of array points, centroids - pointer to an array of centroids, 
           clust_num - number of clusters, dim - number of axes, iterat_num - max number of iterations, 
           indexes - pointer to an array of centroids' indexes, indec - pointer to an array of clusters' numbers for each point
*/
void k_means(point * pnts, int N, point * centroids, int clust_num, int dim, int iterat_num, int * indexes, int * indec);

/*
overview: finds all the neighbours of the point considering the value of epsilon
arguments: pnt - pointer to the root of kd-tree of points, processing - point for which the neigbour points are being found,
           neigh_count - pointer to an int storing the number of neighbours of a point, indexes - pointer to an array of indexes
           of point's neighbours, epsilon - the maximum distance for which the points are still considered as neighbours, 
           ind - index of the dimention used for comparison, dim - number of axes
*/
void all_neighbours(Node * pnt, point processing, int * neigh_count, int * indexes, double epsilon, int ind, int dim);

/*
overview: makes a cluster of points in dbscan algorithm
arguments: root - pointer to the root of kd-tree of points, pnts - pointer to an array of points, N - size of an array pnts,
           neigh_count - pointer to an int storing the number of neighbours of a point, indexes - pointer to an array of indexes
           of point's neighbours, epsilon - the maximum distance for which the points are still considered as neighbours, 
           min_density - minimal density and the minimum number of neighbours for which the point is considered as main, 
           visit - pointer to an array for flags of visits for points, dim - number of axes, indec - pointer to an array of 
           clusters' numbers for each point, clust_num - number of clusters, current_ind - index of a point being processed
*/
void making_cluster(Node * root, point * pnts, int N, int * neigh_count, int * indexes, double epsilon, int min_density, 
    int * visit, int dim, int * indec, int clust_num, int current_ind);

/*
overview: performs a dbscan algorithm dividing pnts into clusters and finding outliers
arguments: pnts - pointer to an array of points, N - size of an array pnts, epsilon - the maximum distance for which the points are 
           still considered as neighbours, min_density - minimal density and the minimum number of neighbours for which the point is 
           considered as main, dim - number of axes, indec - pointer to an array of clusters' numbers for each point
return: the number of clusters
*/
int dbscan(point * pnts, int N, double epsilon, int min_density, int dim, int * indec);

/*
overview: counts the coordinates for each centroid in the dbscan algorithm
arguments: centres - pointer to a pointer to an array of centres, pnts - pointer to an array of points, clust_count - number of clusters + 1,
           N - size of an array pnts, dim - number of axes, indec - pointer to an array of clusters' numbers for each point
*/
void clusters_for_dbscan(point ** centres, point * pnts, int clust_count, int N, int dim, int * indec);

/*
overview: reads the data - first string of coordinates (each coordinate is a string without spaces), and then N points of the dim
arguments: N - number of points to read, dim - pointer to an int storing the number of axes
return: pointer to an array of points
*/
point *  data_set_read(int N, int * dim);

/*
overview: prints the centroids
arguments: centroids - pointer to an array of centroids, clust_num - size of an array centroids, dim - number of axes, flag - 
           flag indicating the algorithm used - 0 - k-means, 1 - dbscan
*/
void printing_centroids(point * centroids, int clust_num, int dim, int flag);

/*
overview: prints all of the points with the number of the cluster they belong to (or also outlier for some points in dbscan)
arguments: points - pointer to an array of points, N - size of an array points, dim - number of axes, indec - pointer to an array of 
           clusters' numbers for each point, flag - flag indicating the algorithm used - 0 - k-means, 1 - dbscan
*/
void printing_points(point * points, int N, int dim, int * indec, int flag);