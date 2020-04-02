//lab1_pthread.cpp

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <algorithm>

using namespace std;

pthread_mutex_t lock;

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

struct thread_data{

	int start;
	int end;
	int K;
	bool* done;
	Cluster** clusters;
	Point** points;

};

void* nearest_cluster(void* threadarg){
    
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
    Point *points;
    points = *my_data->points;
    Cluster *clusters;
    clusters = *my_data->clusters;

	for(int i = my_data->start; i< my_data->end; i++)
    {
        
        int current_cid = points[i].cid;
	    
	    double sum = 0.0, min_dist;
		
		int nearest_cid;

	    sum += pow(clusters[0].cx - (float)points[i].x, 2.0);
	    sum += pow(clusters[0].cy - (float)points[i].y, 2.0);
	    sum += pow(clusters[0].cz - (float)points[i].z, 2.0);
	    
	    min_dist = sqrt(sum);
	    nearest_cid = clusters[0].cid;

	    for(int k = 1; k < my_data->K; k++)
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

        pthread_mutex_lock(&lock);

        if(current_cid != nearest_cid)
        {
            if(current_cid != 0){
                for(int j=0; j<my_data->K; j++){
                    if(clusters[j].cid == current_cid){
                        clusters[j].removePoint(points[i].pid);
                    }
                }
            }

            for(int j=0; j<my_data->K; j++){
                if(clusters[j].cid == nearest_cid){
                    clusters[j].addPoint(points[i]);
                }
            }
            
            points[i].cid = nearest_cid;
            *(my_data->done) = false;
        }
        pthread_mutex_unlock(&lock);

    }

}

void kmeans_pthread(int num_threads,
					int N,
					int K,
					int* data_points,
					int** data_point_cluster,
					float** centroids,
					int* num_iterations
					){

    float *cent= (float*)malloc((K+1)*3*1000*sizeof(float));
    int *dpc = (int*)malloc((N+1)*4*sizeof(int));

	Cluster *clusters;
    Cluster tempc[K];
    clusters = &tempc[0];
	Point *points;
    Point temp[N];
    points = &temp[0];


	pthread_t threads[num_threads];
    struct thread_data td[num_threads];

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
    	cent[3*i]=clusters[i].cx;
    	cent[(3*i)+1]=clusters[i].cy;
    	cent[(3*i)+2]=clusters[i].cz;
    }

    // starting kmeans clustering

    int iteration =0;

    int iter = 0;
    while(iteration<50)
    {
        
        // iteration number - iter

        bool done = true;

        // adding all points to their nearest cluster

        int rc;
        
        for(int yy=0; yy<num_threads; yy++){
        	td[yy] = { .start = (N/num_threads)*yy, .end = (N/num_threads)*(yy+1), .K = K, .done = &done};
        	td[yy].points = &points;
            td[yy].clusters = &clusters;
        }

        for(int x=0; x<num_threads; x++){
        	rc = pthread_create(&threads[x], NULL, &nearest_cluster, (void *) &td[x]);
        	//pthread_join(threads[x],NULL);
        	if(rc)
        		cerr<<"nhi bna";
        }

        for(int x=0; x<num_threads; x++){
        	pthread_join(threads[x],NULL);
        }

        
        // recalculating centroid of each cluster

        for(int i = 0; i < K; i++)
        {
            int ClusterSize = clusters[i].points.size();

            float sumx = 0.0;
            float sumy = 0.0;
            float sumz = 0.0;

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
    	cent[(iter*K)+(3*i)]=clusters[i].cx;
    	cent[(iter*K)+(3*i)+1]=clusters[i].cy;
    	cent[(iter*K)+(3*i)+2]=clusters[i].cz;
    	}

        iteration++;
    }

    for(int i=0; i<N; i++){

		dpc[4*i] = points[i].x;
		dpc[(4*i)+1] = points[i].y;
		dpc[(4*i)+2] = points[i].z;
		dpc[(4*i)+3] = points[i].cid;
	
	}


	*data_point_cluster = dpc;
	*centroids = cent;
	*num_iterations = iter;

}