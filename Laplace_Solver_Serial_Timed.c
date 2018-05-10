#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

// Paramaters
#define X_Min -1.
#define X_Max 1.
#define Y_Min -1.
#define Y_Max 1.
#define PI 3.14159265358979323846
#define N_Mesh 500
#define Change_Threshold (.5/(double)(N_Mesh*N_Mesh))
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
void Update(double *U_Odd, double *U_Even, int Iter);		// Used to make U^(k+1) from U^(k)
double Maximum_Change(double *U_Odd, double *U_Even);	// Calculates average change at meshpoints
void Save_To_File(double *U_k);						// Saves the result in a file.


int main() {
	///////////////////////////////////////////////////////////////////////////
	// Set up timing variables
	clock_t timer, runtime_timer;
	float t_Alloc, t_IC, t_BC, t_Iter, t_Save, t_runtime;

	///////////////////////////////////////////////////////////////////////////
	// Allocation and Initialization: 	

	// Allocate U_Odd and U_Even array. These will store the value of u^(k) for
	// successive iterations. As their names would suggest, U_Odd is used to store
	// Odd iterations while U_Even is used to store even iterations. 
	// In general, if k is odd U^(k+1) is stored in U_Odd and calculated from U_Even
	// (which should be storing U^(k)). Likewise, if k is even then U^(k+1) is 
	// stored in U_Even and calculated from U_Odd (which should be storing U^(k)).
	// Each array is populated with the boundary conditions (since we will be using
	// Both matricies to update). Only U_Even is set up with initial conditions 
	// This is because the first iteration will be for iter = 1, meaning that the
	// 1st iteration is stored in U_Odd and calulated from U_Even.
	double *U_Even, *U_Odd;

	runtime_timer = clock();
	timer = clock();

	U_Odd = (double *)malloc(sizeof(double)*(N_Mesh+2)*(N_Mesh+2));
	U_Even = (double *)malloc(sizeof(double)*(N_Mesh+2)*(N_Mesh+2));

	t_Alloc = (float)((clock() - timer)/((float)CLOCKS_PER_SEC));


		// Set up initial conditions, boundary conditions. Note that we run the
		// Initial conditions on U_Even because Update first moves the interior
		// elements of U_Even to U_Odd (see update function). By contrast, we
		// need BC's on both U_Odd and U_Even because both will be used to 
		// calculate updated iterations

	timer = clock();
	Initial_Conditions(U_Even); 
	t_IC = (float)((clock() - timer)/((float)CLOCKS_PER_SEC));

	timer = clock();
	Boundary_Conditions(U_Even);
	Boundary_Conditions(U_Odd);
	t_BC = (float)((clock() - timer)/((float)CLOCKS_PER_SEC));


	///////////////////////////////////////////////////////////////////////////
	// Update loop:

	// Update matrix until Maximum change falls below threshold
	double Max_Change = 1;
	int Iter = 0;

	timer = clock();

	do {
		// Increment number of iterations 
		Iter++;

		// Perform next iteration
		Update(U_Odd, U_Even, Iter);

		// Find maximum change from this iteration
		Max_Change = Maximum_Change(U_Odd, U_Even);
	} while(Max_Change > Change_Threshold);

	t_Iter = (float)((clock() - timer)/((float)CLOCKS_PER_SEC));

	// Save results (Note, which matrix U^(k+1) is stored in depends on
	// the value of Iter)
	timer = clock();
	// If an Even number of iterations occured, U^(k+1) is in U_Even
	if(Iter%2 == 0) {
		Save_To_File(U_Even);
	} // if(Iter%2 == 0) {

	// Otherwise, U^(k+1) is stored in U_Odd
	else {
		Save_To_File(U_Odd);
	} // else
	t_Save = (float)((clock() - timer)/((float)CLOCKS_PER_SEC));

	t_runtime = ((float)(clock() - runtime_timer))/((float)CLOCKS_PER_SEC);

	// Print timing results
	printf("\t\t -- Paramaters --\n\n");
	printf("Number of meshpoints         ::    %d\n",N_Mesh);
	printf("Change threshold             ::    %f\n",Change_Threshold);
	printf("Number of iterations needed  ::    %d\n",Iter);
	printf("\n\t\t -- Timing data --\n\n");
	printf("Time to alloc U_Odd,U_Even   ::    %.2e (s)\n",t_Alloc);
	printf("Time to set IC's             ::    %.2e (s)\n",t_IC);
	printf("Time to set BC's             ::    %.2e (s)\n",t_BC);
	printf("Time for iterations          ::    %.2e (s)\n",t_Iter);
	printf("Time to save to file         ::    %.2e (s)\n",t_Save);
	printf("Total runtime                ::    %.2e (s)\n",t_runtime);

	return 0;
} // int main()

void Initial_Conditions(double *U_0) {
	// Here, we populate the interior elements of the U_0 matrix.
	// We do not populate the boundary elements, since these are
	// taken care of by the boundary conditions. Finally, recall that
	// U_Even and U_Odd are both of dimension N_Mesh+2. Therefore, the
	// Interior elements have indicies 1,2,...N_Mesh.

	// For simplicity, we will populate the interior elements with a value of 0
	unsigned int i,j;

	for(i = 1; i < N_Rows-1; i++) {
		for(j = 1; j < N_Cols-1; j++) {
			U_0[i*N_Cols + j] = 0;
		} // for(j = 0; j < N_Cols; j++) {
	} // for(i = 0; i < N_Rows; i++) {
} // void Initial_Conditions(double *U_0) {

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
	} // 	for(i = 1; i < N_Mesh+1; i+=) {
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


void Update(double *U_Odd, double *U_Even, int Iter) {
	unsigned int i,j;
	double *U_k, *U_kp1;

	// If Iter is odd, then U_Odd is U^(k+1) and U_Even is U^(k). 
	// Thus, want to build U_Odd from U_Even 
	if( Iter%2 == 1) {
		U_kp1 = U_Odd;
		U_k = U_Even;
	} // 	if( Iter%2 == 1) {

	// If Iter is Even, then U_Even is U^(k+1) and U_Odd is U^(k). 
	// Thus, want to build U_Even from U_Odd
	else {
		U_kp1 = U_Even;
		U_k = U_Odd;
	} // else 

	// Now, calculate the new U^(k+1) values, store them in the interior
	// Elements of U_kp1.
	for(i = 1; i < N_Rows-1; i++) {
		for(j = 1; j < N_Cols-1; j++) {
			U_kp1[i*(N_Cols) + j] = (1./(2.*dx*dx + 2.*dy*dy))*
			                                  ( (dx*dx)*
				                              ( U_k[(i+1)*N_Cols + j    ]   + 
			                                    U_k[(i-1)*N_Cols + j    ])  +
			                                  (dy*dy) * 
			                                  ( U_k[i*N_Cols     + (j-1)]   + 
			                                    U_k[i*N_Cols     + (j+1)]) );
		} // for(j = 1; j < N_Rows-1; j++) {
	} // for(int i = 1; i < N_Cols-1; i++) {
} // void Update(double *U_Odd, double *U_Even, int Iter) {

double Maximum_Change(double *U_Odd, double *U_Even) {
// This finds the maximum change the value of u betweek the kth and 
	// k+1th itteration. This is done by finding the maximum value of 
	// abs(U_Odd(i,j) - U_Even(i,j)) for all interior elements. 
	// It should be noted that we don't need to know which one corresonds
	// to U^(kp1) vs U^(k) since we are finding the absolute value of the
	// difference between each cell of the two matricies.

	double Max_Change = 0;
	double Current_MeshPoint_Change;
	unsigned int i,j;
	for(i = 1; i < N_Rows-1; i++) {
		for(j = 1; j < N_Cols-1; j++) {
			Current_MeshPoint_Change = fabs(U_Odd[i*N_Cols + j] - U_Even[i*N_Cols+j]);

			if (Current_MeshPoint_Change > Max_Change) {
				Max_Change = Current_MeshPoint_Change;
			} // if (Current_Cell_Change > Max_Change) {
		} // for(j = 1; j < N_Cols-1; j++) {
	} // for(i = 1; i < N_Rows-1; i++) {

	return Max_Change;
} // double Avg_Cange(double *U_Odd, double *U_Even) {

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