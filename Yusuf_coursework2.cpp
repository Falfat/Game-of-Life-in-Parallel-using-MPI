#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>

//Defining tag numbers to be used for sends and recieves.
#define Tag_left 0
#define Tag_right 0
#define Tag_top 1
#define Tag_down 1
#define Tag_top_left 3
#define Tag_bottom_left 2
#define Tag_top_right 2
#define Tag_bottom_right 3

//This is to specify if the choice of "periodic" or "non periodic" is preferred
/*
periodic	true : is for a periodic system choice
periodic	false : is for a non periodic system choice
*/
#define periodic	true

using namespace std;

//globally declare variables.
/*
id: is the processor identification number.
p: total number of processors.
p_rows: number of rows in a specific processor
p_cols: number of columns in a specific processor
*/
int id, p;
int p_rows, p_cols;

//declare domain and new_domain
bool **domain;
bool **new_domain;

//globally declare the global imax and jmax.
int cells_imax;
int cells_jmax;
/*
imax_local = total number of cells on the x-axis.
jmax_local = total number of cells on the y-axis.
*/
int imax_local, jmax_local;

//specify maximum number of steps
int max_steps = 100;

/*
div_x: vector to store the local column sizes for the processors
div_y: vector to store the local row sizes for the processors
*/
vector<int>div_x;
vector<int>div_y;

/*
my_row = the row of a processor in the domain
my_column = the column of a processor in the domain
*/
int my_row, my_column;


//function to get dimension of processors for a specified number of processors
//creates domains for processors.
void find_dimensions(int p, int &rows, int &columns)		//A bit brute force - this can definitely be made more efficient!
{
	//let minimum gap equals total number of processors.
	int min_gap = p;

	//loop until the minimum total processors can be divided to : sqrt(p)
	for (int i = 1; i <= sqrt(p) + 1; i++)
	{
		//This loop is executed if the remainder after dividing the total number of processor by iteration
		//number is 0.
		if (p%i == 0)
		{
			//calculate gap
			int gap = abs(p / i - i);

			//executed if the calculated gap is less than minimum gap
			if (gap < min_gap)
			{
				//let minimum gap equals gap
				min_gap = gap;
				//number of rows equals iteration number
				rows = i;
				//number of colums equal total number of processor divided by iteration number
				columns = p / i;
			}
		}
	}
	//print out the dimensions of the domains we have created.
	if (id == 0)
		cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
}


//function to determine the cell size(x,y) in the local domain of a processor
void cells_for_grids(int rows, int columns, int imax, int jmax)
{
	//give sizes of values vectors should store
	div_x.resize(columns);
	div_y.resize(rows);

	//declare remaining cells left in both x and y directions/
	int remaining_y = imax;
	int remaining_x = jmax;
	//loop to allocate - in the x-direction- cells to processors
	for (int i = 0; i < rows; i++)
	{
		//allocate for each processor down the grid
		int allocated = remaining_y / (rows - i);
		//store the allocated cells for the processor in the vector
		div_y[i] = allocated;
		//deduct the allocated from the global imax
		remaining_y -= allocated;
	}

	//loop to allocate - in the y-direction- cells to processors
	for (int j = 0; j < columns; j++)
	{
		//allocate for each processor across the grid
		int allocated = remaining_x / (columns - j);
		//store the allocated cells for the processor in the vector
		div_x[j] = allocated;
		//deduct the allocated from the global jmax 
		remaining_x -= allocated;
	}

}

//function to convert processors id to it's index on the domain
/*
for instance processor 0's index will be (0,0)
*/
void id_to_index(int id, int &id_row, int &id_column)
{
	//get the column index for the processor id as specified in the argument
	id_column = id % p_cols;
	//get the row idex for the processor id as specified in the argument
	id_row = id / p_cols;
}

//function to determine processor's id from its index.
/*
For instance, index (0,0) will give id = 0.
*/
int id_from_index(int id_row, int id_column)
{
	//if periodic is true
	if (periodic)
	{
		//return processor's id for a periodic domain
		return ((id_row + p_rows) % p_rows) * p_cols + (id_column + p_cols) % p_cols;
	}
	//executed for a non-periodic domain.
	else
	{
		//check if the processor's row id is out of bounds.
		if (id_row >= p_rows || id_row < 0)
			return -1;
		//check if the processor's column id is out of bounds.
		if (id_column >= p_cols || id_column < 0)
			return -1;
		//return processor's id
		return id_row * p_cols + id_column;
	}
}

//function to setup peer to peer communication among processors.
void comms(void)
{
	//declare count and initialize it to be 0
	int cnt = 0;

	//declare and allocate memory for sending and recieving.
	bool *send_right = new bool[imax_local];
	bool *send_left = new bool[imax_local];
	bool *recv_left = new bool[imax_local];
	bool *recv_right = new bool[imax_local];

	//fill the values to be sent left.
	for (int i = 0; i < imax_local; i++)
		send_left[i] = domain[i + 1][1];


	//fill the values to be sent right.
	for (int i = 0; i < imax_local; i++)
		send_right[i] = domain[i + 1][jmax_local];


	//declare and allocate memory for request to send and recieve
	MPI_Request *request = new MPI_Request[8 * 2];


	//send data up
	int com_id = id_from_index(my_row - 1, my_column);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[1][1], jmax_local, MPI_BYTE, com_id, Tag_top, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[0][1], jmax_local, MPI_BYTE, com_id, Tag_down, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process(this can occur when using a few processes)
	else if (com_id == id)
	{
		//copy data back to same process' padding down
		for (int i = 1; i <= jmax_local; i++)
			domain[imax_local + 1][i] = domain[1][i];
	}


	//send data down
	com_id = id_from_index(my_row + 1, my_column);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[imax_local][1], jmax_local, MPI_BYTE, com_id, Tag_down, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[imax_local + 1][1], jmax_local, MPI_BYTE, com_id, Tag_top, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding up
		for (int i = 1; i <= jmax_local; i++)
			domain[0][i] = domain[imax_local][i];
	}

	//send data to the right
	com_id = id_from_index(my_row, my_column + 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(send_right, imax_local, MPI_BYTE, com_id, Tag_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(recv_right, imax_local, MPI_BYTE, com_id, Tag_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding left
		for (int i = 1; i < imax_local; i++)
			domain[i][jmax_local] = domain[i][0];
	}

	//send data to the left
	com_id = id_from_index(my_row, my_column - 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(send_left, imax_local, MPI_BYTE, com_id, Tag_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(recv_left, imax_local, MPI_BYTE, com_id, Tag_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding right
		for (int i = 1; i < imax_local; i++)
			domain[i][1] = domain[i][jmax_local + 1];
	}

	//send data top-right
	com_id = id_from_index(my_row - 1, my_column + 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[1][jmax_local], 1, MPI_BYTE, com_id, Tag_top_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[imax_local + 1][0], 1, MPI_BYTE, com_id, Tag_bottom_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding bottom-left
		domain[imax_local + 1][0] = domain[1][jmax_local];
	}

	//send data top-left
	com_id = id_from_index(my_row - 1, my_column - 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[1][1], 1, MPI_BYTE, com_id, Tag_top_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[imax_local + 1][jmax_local + 1], 1, MPI_BYTE, com_id, Tag_bottom_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding bottom-right
		domain[imax_local + 1][jmax_local + 1] = domain[1][1];
	}


	//send data down-left
	com_id = id_from_index(my_row + 1, my_column - 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[imax_local][1], 1, MPI_BYTE, com_id, Tag_bottom_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[0][jmax_local + 1], 1, MPI_BYTE, com_id, Tag_top_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding top-right
		domain[0][jmax_local + 1] = domain[imax_local][1];
	}


	//send data down-right
	com_id = id_from_index(my_row + 1, my_column + 1);

	//send to nighbour up and not itself. if not periodic, dont send out of bounds
	if (com_id >= 0 && com_id != id)
	{
		MPI_Isend(&domain[imax_local][jmax_local], 1, MPI_BYTE, com_id, Tag_bottom_right, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
		MPI_Irecv(&domain[0][0], 1, MPI_BYTE, com_id, Tag_top_left, MPI_COMM_WORLD, &request[cnt]);
		cnt++;
	}
	//if the matching periodic boundary is on the same process
	else if (com_id == id)
	{
		//copy data back to same process' padding top-left
		domain[0][0] = domain[imax_local][jmax_local];
	}

	//wait for all sends and recieves to complete
	MPI_Waitall(cnt, request, MPI_STATUSES_IGNORE);

	
	//copy into padding after recieving to the left
	for (int i = 0; i < imax_local; i++)
		domain[i + 1][0] = recv_left[i];

	//copy into padding after recieving to the right
	for (int i = 0; i < imax_local; i++)
		domain[i + 1][jmax_local + 1] = recv_right[i];


	//delete pointers
	delete[] request;
	delete[] send_left;
	delete[] send_right;
	delete[] recv_left;
	delete[] recv_right;

}

//function to check number of neighbours alive
int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++) //loop around all the eight neighbours
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0) //avoid counting yourself
			{
				ix = i + ii; //specify the row to be checked
				jx = j + jj; //specify the column to be checked
				//if alive
				if (domain[ix][jx]) cnt++; //increase count
			}
	return cnt; //return total number of neighbours alive.
}

//function to do iteration
void do_iteration(void)
{
	//loop through the local domain of the processor
	for (int i = 1; i <= imax_local; i++)
		for (int j = 1; j <= jmax_local; j++)
		{
			//let the new status (dead/alive) of the cell equal previous status
			new_domain[i][j] = domain[i][j];
			//get the number of neighbours of the cell alive
			int num_n = num_neighbours(i, j);
			if (domain[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_domain[i][j] = false;
			}
			//if the cell has fewer than 2 neighbours and 4 or more neighbours
			else if (num_n == 3) new_domain[i][j] = true;
		}
	//update the domain before the next iteration.
	bool **temp;
	temp = domain;
	domain = new_domain;
	new_domain = temp;
}


//function to print the results from each processor after every iteration.
void grid_to_file(int it)
{
	stringstream fname;
	fstream f1;
	fname << "output" << "_" << id << "_" << it << ".txt";
	//open file
	f1.open(fname.str().c_str(), ios_base::out);
	f1 << "imax local is" << endl;
	f1 << imax_local << endl;
	f1 << "jmax is local" << endl;
	f1 << jmax_local << endl;
	f1 << "my rows is" << endl;
	f1 << my_row << endl;
	f1 << "my cols is" << endl;
	f1 << my_column << endl;
	f1 << "number of rows " << endl;
	f1 << p_rows << endl;
	f1 << "number of cols " << endl;
	f1 << p_cols << endl;
	f1 << "number of cores " << endl;
	f1 << p << endl;
	f1 << "global imax is " << endl;
	f1 << cells_imax << endl;
	f1 << "global jmax is " << endl;
	f1 << cells_jmax << endl;
	//loop through the processor's local domain and print without padding.
	for (int i = 1; i <= imax_local; i++)
	{
		for (int j = 1; j <= jmax_local; j++)
			f1 << domain[i][j] << "\t";
		f1 << endl;
	}

	//close file
	f1.close();
}


int main(int argc, char *argv[])
{
	//setup MPI communications
	MPI_Init(&argc, &argv);
	//Reads the rank (number) of the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	//Reads the total number of processes that have been assigned
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//seeding for generating random values
	srand(time(NULL) + id * 10);

	//initialize the size of the global grid
	cells_imax = 100;
	cells_jmax = 100;
	

	//find the dimensions the total number of processors will be divided into.
	find_dimensions(p, p_rows, p_cols);

	//determine the cell size(x, y) in the local domain of each processor
	cells_for_grids(p_rows, p_cols, cells_imax, cells_jmax);

	//create the indices of the processors in the global domain from their id
	id_to_index(id, my_row, my_column);

	//each processor knows its  local imax and jmax
	imax_local = div_y[my_row];
	jmax_local = div_x[my_column];

	//allocate memory for the local domain i_max of each processor including buffers
	domain = new bool*[imax_local + 2];
	//allocate memory for the local new_domain i_max of each processor including buffers. This is needed for iteration.
	new_domain = new bool*[imax_local + 2];

	//allocate memory for the local domain j_max of each processor including buffers
	for (int i = 0; i < imax_local + 2; i++)
		domain[i] = new bool[jmax_local + 2];

	//allocate memory for the local new_domain j_max of each processor including buffers. This is needed for iteration.
	for (int i = 0; i < imax_local + 2; i++)
		new_domain[i] = new bool[jmax_local + 2];

	
	//initialise domain with random values of 0 and 1 only. 1 means alive and 0 means dead.
	//loop round the loacal domain of each processor, excluding the buffer.
	for (int i = 1; i < imax_local + 1; i++)
		for (int j = 1; j < jmax_local + 1; j++)
		{
			domain[i][j] = (rand() % 2);
			new_domain[i][j] = 0;
		}

	for (int n = 0; n < max_steps; n++)
	{
		//set up communication with neighbouring processors (peer-peer communication)
		comms();
		//print out the number of iterations concluded from first processor
		if (id == 0)
			cout << "it: " << n << endl;
		do_iteration(); //do iterations to implement the game of live rules.

		grid_to_file(n); //output file for each iteration

	}
	
	//delete pointers to pointers for domain
	for (int i = 0; i < imax_local + 2; i++)
		delete[] domain[i];
	delete[] domain;

	//delete pointers to pointers for new_domain
	for (int i = 0; i < imax_local + 2; i++)
		delete[] new_domain[i];
	delete[] new_domain;

	MPI_Finalize();

	//system("pause");
}