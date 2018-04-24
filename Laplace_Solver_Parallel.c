#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "omp.h"

// Paramaters
#define X_Min -1.
#define X_Max 1.
#define Y_Min -1.
#define Y_Max 1.
#define PI 3.14159265358979323846
#define N_Mesh 10
#define Change_Threshold (.001/(double)N_Mesh)
const unsigned int N_Rows = N_Mesh+2;
const unsigned int N_Cols = N_Mesh+2;

// Determine dx,dy (spacing between successive gridpoints in the x and y directions
// Note that if there are D_Div grid points in a given direction, then there are
// N_Mesh-1 jumps (of size dx) that have to be made to get from X_Min to X_Max
#define dx ((X_Max - X_Min)/((double)N_Mesh-1.))
#define dy ((Y_Max - Y_Min)/((double)N_Mesh-1.))

// prototypes
void Initial_Conditions(double *U_0);				// Used to set IC's
void Boundary_Conditions(double *U_0);				// Used to set BC's
	double Left_BC(double y);						// BC function for left boundary
	double Right_BC(double y);						// BC function for right boundary
	double Top_BC(double y);						// BC function for top boundary
	double Bottom_BC(double y);						// BC function for bottom boundary
void Update(double *U_k, double *U_kp1);			// Used to make U^(k+1) from U^(k)
void Maximum_Change(double *U_k, double *U_kp1, double *Max_Change);	// Calculates average change at meshpoints
void Save_To_File(double *U_k);						// Saves the result in a file.
	
int main() {
	///////////////////////////////////////////////////////////////////////////
	// Set up variables
	double timer, t_Alloc, t_IC, t_BC, t_Iter, t_Save;

	double Max_Change = 1;
	int iterations = 0;

	///////////////////////////////////////////////////////////////////////////
	// Allocation and initlaization: 	

	// Allocate U_k and U_kp1 array. U_k will store the values for u^(k). 
	// During the update cycle this array will be read only. U_kp1 stores the 
	// values for U^(k+1). We make U_k of dimension (N_Mesh+2)x(N_Mesh+2)
	// So that we can include boundary elements. We also make U_kp1 the same size
	// to make access more intuitive. 
	double *U_k, *U_kp1;

	timer = omp_get_wtime();

	U_k = (double *)malloc(sizeof(double)*(N_Mesh+2)*(N_Mesh+2));
	U_kp1 = (double *)malloc(sizeof(double)*(N_Mesh+2)*(N_Mesh+2));

	t_Alloc = omp_get_wtime() - timer;

	// Parallel variables
	int n_procs = omp_get_num_procs();
	omp_set_num_threads(n_procs;

	#pragma omp parallel default(none) firstprivate(n_procs, N_Cols, N_Rows) shared(iterations, Max_Change, U_k, U_kp1,t_Alloc, t_IC, t_BC, t_Iter, t_Save, timer)
	{

		// Set up initial conditions, boundary conditions. Note that we run the
		// Initial conditions on U_kp1 because Update first moves the interior
		// elements of U_kp1 to U_k (see update function)

		#pragma omp master
		{ timer = omp_get_wtime(); }
		Initial_Conditions(U_kp1); 
		#pragma omp master
		{ t_IC = omp_get_wtime() - timer; }

		#pragma omp master
		{ timer = omp_get_wtime(); }
		Boundary_Conditions(U_k);
		#pragma omp master
		{ t_BC = omp_get_wtime() - timer; }

		///////////////////////////////////////////////////////////////////////////
		// Update loop:

		// Update matrix until Maximum change falls below threshold
		#pragma omp master
		{ timer = omp_get_wtime(); }

		while(Max_Change > Change_Threshold) {
			// Perform next iteration
			Update(U_k, U_kp1);

			// Find maximum change from this iteration
			Maximum_Change(U_k, U_kp1, &Max_Change);

			// Incremeent number of iterations
			#pragma omp master
			{ iterations++; }
		} // while(Max_Change > Change_Threshold) {
		
		// Update once more, We do this so that u_kp1 from the final iteration
		// is moved to the U_k matrix.
		Update(U_k, U_kp1);

		#pragma omp master
		{ t_Iter = omp_get_wtime() - timer; }

		// Save results. We only want one copy (so just the master thread does it.)
		#pragma omp master 
		{
			timer = omp_get_wtime();
			Save_To_File(U_k);
			t_Save = omp_get_wtime() - timer;
		} // #pragma omp master 

	} // #pragma omp parallel default(none) firstprivate(n_procs, N_Cols, N_Rows) shared(iterations, Max_Change, U_k, U_kp1,t_Alloc, t_IC, t_BC, t_Iter, t_Save, timer)

	// Print timing results
	printf("\t\t -- Paramaters --\n\n");
	printf("Number of procs              ::\t %d\n",omp_get_num_procs());
	printf("Number of threads            ::\t %d\n",omp_get_max_threads());
	printf("Number of meshpoints         ::\t %d\n",N_Mesh);
	printf("Change threshold             ::\t %f\n",Change_Threshold);
	printf("Number of iterations needed  ::\t %d\n",iterations);
	printf("\n\t\t -- Timing data --\n\n");
	printf("Time to alloc U_k,U_kp1      ::\t %.2e (s)\n",t_Alloc);
	printf("Time to set IC's             ::\t %.2e (s)\n",t_IC);
	printf("Time to set BC's             ::\t %.2e (s)\n",t_BC);
	printf("Time for iterations          ::\t %.2e (s)\n",t_Iter);
	printf("Time to save to file         ::\t %.2e (s)\n",t_Save);

	return 0;
} // int main()

void Initial_Conditions(double *U_0) {
	// Here, we populate the interior elements of the U_0 matrix.
	// We do not populate the boundary elements, since these are
	// taken care of by the boundary conditions. Finally, recall that
	// U_K and U_kp1 are both of dimension N_Mesh+2. Therefore, the
	// Interior elements have indicies 1,2,...N_Mesh.

	// For simplicity, we will populate the interior elements with a value of 0
	unsigned int i,j;

	#pragma omp for schedule(static) private(i,j)
		for(i = 1; i < N_Rows-1; i++) {
			for(j = 1; j < N_Cols-1; j++) {
				U_0[i*N_Cols + j] = 0;
			} // for(j = 0; j < N_Cols; j++) {
		} // for(i = 0; i < N_Rows; i++) {
} // void Initial_Conditions(double *U_k) {

void Boundary_Conditions(double *U_0) {
	// Here we populate the boundary elements of U_0. 
	// This function can make use of X_Max, X_Min, Y_Max, and Y_Min you
	// want to use a function to define the boundary.

	// It should be noted that we do not need to populate the corner elements
	// Since these will not be used to update any interior elements. 

	// Here we are populating the bottom edge of our matrix with a value of 1
	// The rest of the boundary is given a value of 0

	// Set up x, y
	double x,y;

	unsigned int i,j;

	#pragma omp for schedule(static) private(i)
		for(i = 1; i < N_Mesh+1; i++) {
			// Set x, y based on position
			x = (i-1)*dx + X_Min;
			y = Y_Max - (i-1)*dy;		// Note it is easier to decrease y as we move down the rows 

			//////////////////////////////////////////////////////////////////////////////
			// Vertical boundaries:

			// Left Edge:
			// This is simply the first column of U_0. Each column is separated by
			// N_Col elements in the U_0 array.
			U_0[(N_Cols)*i] = Left_BC(y);

			// Right Edge:
			// This is simply the last column of U_0. The Last column occurs at 
			// indicies (N_Cols-1), (N_Cols-1)+1*(N_Cols),... (N_Cols-1)+(N_Rows-1)*(N_Cols)
			U_0[(N_Cols-1)+(N_Cols)*i] = Right_BC(y);

			//////////////////////////////////////////////////////////////////////////////
			// Vertical boundaries:

			// Top Edge:
			// This is simply the first row of A.
			U_0[i] = Top_BC(x);

			// Bottom edge:
			// This is the bottom row of U_0. To get there, we eed to skip the first
			//  N_Mesh+1 rows each row has has N_Cols columns
			U_0[(N_Rows-1)*(N_Cols) + i] = Bottom_BC(x);
		} // for(i = 1; i < N_Mesh+1; i+=) {
} // void Boundary_Conditions(double *U_0) {

double Left_BC(double y) {
	// This function us used to set the left boundary condition function.
	// Since the left boundary is a vertical boundary, this is a function of y
	double LeftBC_y = (1./2.)*(1.-sin(y*2.*PI/(Y_Max-Y_Min)));
	return LeftBC_y;
} // double Left_BC(double y) {

double Right_BC(double y) {
	// This function us used to set the right boundary condition function.
	// Since the right boundary is a vertical boundary, this is a function of y
	double RightBC_y  = (1./2.)*(1.+sin(y*2.*PI/(Y_Max-Y_Min)));
	return RightBC_y;
} // double Right_BC(double y) {

double Top_BC(double x) {
	// This function us used to set the bottom boundary condition function.
	// Since the bottom boundary is a horizontal boundary, this is a function of x
	double TopBC_x = (1./2.)*(1.-sin(x*2.*PI/(X_Max-X_Min)));
	return TopBC_x;
} // double Top_BC(double y) {

double Bottom_BC(double x) {
	// This function us used to set the bottom boundary condition function.
	// Since the bottom boundary is a horizontal boundary, this is a function of x
	double BottomBC_x = (1./2.)*(1.+sin(x*2.*PI/(X_Max-X_Min)));
	return BottomBC_x;
} // double Bottom_BC(double y) {


void Update(double *U_k, double *U_kp1) {
	// First, we need to move the elements of U_kp1 over to U_k (since we have begin
	// the next step)
	unsigned int i,j;

	#pragma omp for schedule(static) private(i,j)
		for(i = 1; i < N_Rows-1; i++) {
			for(j = 1; j < N_Cols-1; j++) {
				U_k[i*(N_Cols) + j] = U_kp1[i*(N_Cols) + j];
			} // for(j = 1; j < N_Mesh+1; j++) {
		} // for(int i = 1; i < N_Mesh+1; i++) {

	// Now, calculate the new U^(k+1) values, store them in the interior
	// Elements of U_kp1.
	#pragma omp for schedule(static) private(i,j)
		for(i = 1; i < N_Rows-1; i++) {
			for(j = 1; j < N_Cols-1; j++) {
				U_kp1[i*(N_Cols) + j] = (1./(2.*dx*dx + 2.*dy*dy))*
			                                  ( (dx*dx)*
				                              ( U_k[(i+1)*N_Cols + j    ]   + 
			                                    U_k[(i-1)*N_Cols + j    ])  +
			                                  (dy*dy) * 
			                                  ( U_k[i*N_Cols     + (j-1)]   + 
			                                    U_k[i*N_Cols     + (j+1)]) );
			} // for(j = 1; j < N_Mesh+1; j++) {
		} // for(int i = 1; i < N_Mesh+1; i++) {

} // void Update(double *U_k, double *U_kp1) {

void Maximum_Change(double *U_k, double *U_kp1, double *Max_Change) {
	// This finds the maximum change the value of u betweek the kth and 
	// k+1th itteration. This is done by finding the maximum value of 
	// abs(U_k(i,j) - U_kp1(i,j)) for all interior elements. 

	double Local_Max_Change = 0;
	*Max_Change = 0;
	double Current_MeshPoint_Change;
	unsigned int i,j;

	// Note, I choose to basically write my own reduction function here. Note that
	// All variables above should be default private since this is an orphaned 
	// worksharing dircetive. 
	#pragma omp for schedule(static) private(i,j)
		for(i = 1; i < N_Rows-1; i++) {
			for(j = 1; j < N_Cols-1; j++) {
				Current_MeshPoint_Change = fabs(U_k[i*N_Cols + j] - U_kp1[i*N_Cols+j]);

				if (Current_MeshPoint_Change > Local_Max_Change) {
					Local_Max_Change = Current_MeshPoint_Change;
				} // if (Current_Cell_Change > Local_Max_Change) {
			} // for(j = 1; j < N_Cols-1; j++) {
		} // for(i = 1; i < N_Rows-1; i++) {

	// Combine the local results. This has the effect of a reduction.
	#pragma omp critical
	{
		if(Local_Max_Change > *Max_Change) {
			*Max_Change = Local_Max_Change;
		} // if(Local_Max_Change > Max_Change) {
	} // #pragma omp critical
	#pragma omp barrier

} // void Avg_Cange(double *U_k, double *U_kp1, double *Max_Change) {

void Save_To_File(double *U_k) {
	// This function saves the results to a matrix.

	// First, open our new file
	FILE* Results = fopen("Results.txt","w");

	unsigned int i,j;
	for(i = 0; i < N_Rows; i++) {
		for(j = 0; j < N_Cols; j++) {
			fprintf(Results,"%.4f",U_k[i*N_Cols + j]);
			
			// If we're on the last colon, we want to print a semicolon and a new line.
			// otherwise we want to print a comma and space
			if(j == N_Cols-1) { fprintf(Results,"\n"); }
			else { fprintf(Results, ", "); }
		}
	}
	fclose(Results);
} // void Print_Matrix(double *Mat, unsigned int num_rows, unsigned int num_cols) {