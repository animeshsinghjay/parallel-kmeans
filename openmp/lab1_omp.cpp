//lab1_omp.cpp

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>

using namespace std;


struct Point{

	int pid; // Point ID 
	int cid; // Cluster ID
	
	int x;	// x coordinate
   	int y;	// y coordinate
	int z;	// z coordinate

};


struct Cluster{

	int cid; // Cluster ID
	
	float cx;	// cluster x coordinate
   	float cy;	// cluster y coordinate
	float cz;	// cluster z coordinate

	vector<Point> points;

	void addPoint(Point p){
    	p.cid = this->cid ;
    	points.push_back(p);
	}

	bool removePoint(int pid){

    	int size = points.size();

    	for(int i = 0; i < size; i++)
        {
            if(points[i].pid == pid)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
	}

};

void kmeans_omp(int num_threads,
                int N,
                int K,
                int* data_points,
                int** data_point_cluster,
                float** centroids,
                int* num_iterations
                ){

    float *cent= (float*)malloc((K+1)*3*1000*sizeof(float));
    int *dpc = (int*)malloc((N+1)*4*sizeof(int));

    omp_set_num_threads(num_threads);

	Cluster clusters[K];
	Point points[N];

	for(int i=0; i<N; i++){

		points[i].pid = i;
		points[i].cid = 0;
		points[i].x = data_points[3*i];
		points[i].y = data_points[(3*i)+1];
		points[i].z = data_points[(3*i)+2];
	
	}

	// initializing clusters

    vector<int> used_pids;

    for(int i=0; i<K; i++)
    {
        while(true)
        {
            int index = rand() % N;

            if(find(used_pids.begin(), used_pids.end(), index) == used_pids.end())
            {
                used_pids.push_back(index);
                points[index].cid = i;
                Cluster cc = { .cid = i , .cx = (float)points[index].x , .cy = (float)points[index].y , .cz = (float)points[index].z };
                cc.addPoint(points[index]);
                clusters[i]=cc;
                break;
            }
        }
    }

    // writing centroids initially

    for(int i=0; i<K; i++){
    	*(cent + 3*i)=clusters[i].cx;
    	*(cent + (3*i)+1)=clusters[i].cy;
    	*(cent + (3*i)+2)=clusters[i].cz;
    }

    // starting kmeans clustering

    int iter = 0;
    while(true)
    {
        
        // iteration number - iter

        bool done = true;

        // adding all points to their nearest cluster
       #pragma omp parallel for num_threads(num_threads)
        for(int i=0; i<N; i++)
        {
            
            int current_cid = points[i].cid;

            double sum = 0.0, min_dist;
            int nearest_cid;

            sum += pow(clusters[0].cx - (float)points[i].x, 2.0);
            sum += pow(clusters[0].cy - (float)points[i].y, 2.0);
            sum += pow(clusters[0].cz - (float)points[i].z, 2.0);

            min_dist = sqrt(sum);
            nearest_cid = clusters[0].cid;

            for(int k = 1; k < K; k++)
            {
                double dist;
                sum = 0.0;

                sum += pow(clusters[k].cx - (float)points[i].x, 2.0);
                sum += pow(clusters[k].cy - (float)points[i].y, 2.0);
                sum += pow(clusters[k].cz - (float)points[i].z, 2.0);

                dist = sqrt(sum);

                if(dist < min_dist)
                {
                    min_dist = dist;
                    nearest_cid = clusters[k].cid;
                }
            }

            #pragma omp critical
            {
            
            if(current_cid != nearest_cid)
            {
                if(current_cid != 0){
                    for(int j=0; j<K; j++){
                        if(clusters[j].cid == current_cid){
                            clusters[j].removePoint(points[i].pid);
                        }
                    }
                }

                for(int j=0; j<K; j++){
                    if(clusters[j].cid == nearest_cid){
                        clusters[j].addPoint(points[i]);
                    }
                }
                
                points[i].cid = nearest_cid;
                done = false;
            }
            }

        }

        // recalculating centroid of each cluster

        for(int i = 0; i < K; i++)
        {
            int ClusterSize = clusters[i].points.size();

            double sumx = 0.0;
            double sumy = 0.0;
            double sumz = 0.0;

            if(ClusterSize > 0)
            {
                for(int p = 0; p < ClusterSize; p++){

                    sumx += clusters[i].points[p].x;
                    sumy += clusters[i].points[p].y;
                    sumz += clusters[i].points[p].z;
                }
                clusters[i].cx = sumx / ClusterSize;
                clusters[i].cy = sumy / ClusterSize;
                clusters[i].cz = sumz / ClusterSize;
            }
        }

        if(done)
        {
            // clustering completed in iteration number - iter
            break;
        }

        iter++;

        // write cluster centers

        for(int i=0; i<K; i++){
    	*(cent + (iter*K)+(3*i))=clusters[i].cx;
    	*(cent + (iter*K)+(3*i)+1)=clusters[i].cy;
    	*(cent + (iter*K)+(3*i)+2)=clusters[i].cz;
    	}

    
    }

    for(int i=0; i<N; i++){

		*(dpc + (4*i)) = points[i].x;
		*(dpc + (4*i)+1) = points[i].y;
		*(dpc + (4*i)+2) = points[i].z;
		*(dpc + (4*i)+3) = points[i].cid;
	
	}


	*data_point_cluster = dpc;
	*centroids = cent;
	*num_iterations = iter;

}