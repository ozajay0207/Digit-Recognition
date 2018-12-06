// test_model_digit.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "fstream"
#include "string"
#include "math.h"
#include "iomanip"

#define CODE_BOOK_SIZE 32
#define FRAME_SIZE 320
#define FRAMES 60
#define N 5
#define M 32
#define T 60
#define DIGITS 9
#define UTTERENCES 30

using namespace std;

ifstream in, in1;
ofstream out, out1;

char* code_book_file = "codebook.txt";

char obs_sequence_file[100] = "observation_sequence_testing.txt" ;

const int total_window = (FRAME_SIZE + (FRAMES - 1) * 80);
long double code_book[CODE_BOOK_SIZE][12] = { 0 };
long double c_prime_arr[T][12] = { 0 };

int frame_cutting_index = 1;
int start_marker = 0;
int end_marker = 0;

int output_model_name = 0;

//Used Files
//char input_file[100] = "Input/184101035_1_24.txt";
char custom_input_file[100] ;
char input_file[100] = "input.txt";

char* normalized_file = "Normalized.txt";
char* silence_file = "silence_file.txt";
char* trimmed_file = "trim.txt";
char* ri_file = "ri_file.txt";
char* ai_file = "ai_file.txt";
char* ci_file = "ci_file.txt";
char* c_prime_file = "c_prime.txt";
char* hamming_file = "Hamming_window.txt";

char a_average_file[100] = "a_i_j_final.txt";
char b_average_file[100] = "b_i_j_final.txt";
char pi_average_file[100] = "pi_final.txt";

char* alpha_file = "alpha.txt";
long double Pobs_model = 0;
long double max_Pobs_model = 0;

long double tokhura_weight[12] = { 1, 3, 7, 13, 19, 22, 25, 33, 42, 50, 56, 61 };
long double tokhura_dist[5] = { 0 };
const int p = 12;
int flag = 0;
int max_sample_index = 0;
int remove_header_counter = 5;
long double hamming_window[320] = { 0 };
long double current_value;
long double sum_samples = 0;
double dc_shift_value = 0;
double normalization_ratio = 1;
long int zcr_count = 0;
long int no_of_samples = 0;
long int no_of_samples1 = 0;
long double max_sample_value = 1;
long double sample_array[total_window] = { 0 };
long double r[12] = { 0 }, k[12] = { 0 }, alpha[13][13] = { 0 }, E[13] = { 0 }, a[13] = { 0 }, c[12] = { 0 }, c_prime[12] = { 0 }, w[12] = { 0 };
int obs[T + 1];

long double a_i_j[N + 1][N + 1];
long double b_i_j[N + 1][M + 1];
long double pi[N + 1];

int count_check = 0;

long double alpha_arr[N+1][T+1];


//To Remove header from the input files
void remove_header(ifstream& in){
	if (flag){
		remove_header_counter = 0;
		string temp = "";
		while (getline(in, temp) && remove_header_counter){
			remove_header_counter--;
		}
		flag = 0;
	}
	else{
		remove_header_counter = 0;
		string temp = "";
		while (getline(in, temp) && remove_header_counter){
			remove_header_counter--;
		}
	}
}

//To get Hamming Windows values to array from the pre-calculated File
void get_hamming_window(){
	int index_count = 0;
	in.open(hamming_file);
	string temp = "";
	while (in >> temp){
		hamming_window[index_count++] = stod(temp);
	}
	in.close();
}

//To calculate DC Shift
void calculate_dc_shift(){
	int count = 0;
	sum_samples = 0;
	no_of_samples = 0;
	in.open(silence_file);
	flag = 1;
	//cout << "\n................Calculating DC shift..................." << endl;
	remove_header(in);
	string temp = "";

	while (in >> temp && !remove_header_counter){
		current_value = stod(temp);
		count++;
		if (count >= 200){
			sum_samples += current_value;
			no_of_samples++;
		}
	}

	//cout << "\nNo of samples:" << no_of_samples << endl;
	//cout << "Sum:" << sum_samples << endl;
	dc_shift_value = sum_samples / no_of_samples;
	//cout << "DC Shift value:" << dc_shift_value << endl;
	in.close();
}

//To calculate normalization ratio
void calculate_normalization_ratio(){
	int index_count = 0;
	max_sample_value = 0;
	max_sample_index = 0;
	no_of_samples1 = 0;
	//cout << "\n............Reading from " << input_file << "................" << endl;
	in.open(input_file);
	//remove_header(in);
	string temp = "";
	//cout << "\n............Calculating Normalization Ratio................" << endl;
	while (in >> temp && !remove_header_counter){
		index_count++;
		no_of_samples1++;
		current_value = stold(temp);

		//Saving maximum index value and maximum sample value
		if (abs(current_value) > max_sample_value){
			max_sample_value = abs(current_value);
			max_sample_index = index_count;
		}
	}
	//cout << "Max Sample value:" << max_sample_value << endl;
	//cout << "Max Sample index:" << max_sample_index << endl;

	//Calculating normalization ratio
	normalization_ratio = 5000.0 / max_sample_value;
	//cout << "Normalization ratio:" << normalization_ratio << endl;
	in.close();
}

//To remove DC Shift and normalize
void dc_normalize(){

	//Initializing Markers as + or - 640 from the max sample value
	start_marker = max_sample_index - (total_window / 2);
	end_marker = max_sample_index + (total_window / 2);

	if (start_marker < 0){
		//cout << "....................................................................................Adjusting Beg........." << endl;
		start_marker = 1;
		end_marker = total_window + 1;
	}

	if (end_marker>no_of_samples1){
		//cout << "....................................................................................Adjusting End........." << endl;
		end_marker = no_of_samples1 + 1;
		start_marker = no_of_samples1 - total_window;
	}

	//cout << "Start Marker:" << start_marker << endl;
	//cout << "End Marker:" << end_marker << endl;


	//cout << "Total Window:" << total_window << endl;

	int index_count = 0;
	int arr_index = 0;
	in.open(input_file);
	remove_header(in);
	string temp = "";
	//cout << "\n........Removing DC shift and Normalizing File.........." << endl;
	out.open(normalized_file);
	out1.open("sample_array.txt");
	//Subtracting DC shift and Multiplying Normalization ratio
	while (in >> temp){
		index_count++;
		current_value = stod(temp);
		current_value = current_value - dc_shift_value;
		current_value = current_value * normalization_ratio;

		if (index_count >= start_marker && index_count < end_marker){
			sample_array[arr_index++] = current_value;
			out1 << sample_array[arr_index - 1] << endl;
		}

		//Writing the Normalized values to file "normalized.txt"
		out << to_string(current_value) << endl;
	}
	//cout << "First:" << sample_array[0] << endl;
	//cout << "Last:" << sample_array[arr_index - 1] << endl;
	//cout << "arr_index" << arr_index << endl;
	out.close();
	out1.close();
	in.close();
}

//Calculating Ri's
void calculate_Ris(){
	//cout << "\n........Writing Ri's to file.........." << endl;
	int count = 0;
	long double first_value = 0;
	long double second_value = 0;
	string temp;

	out.open(ri_file);

	//Calculating Ri and also multiplying with hamming window appropriately
	for (int k = 0; k < FRAMES; k++){
		for (int j = 0; j <= p; j++){
			r[j] = 0;
			for (int i = 0; i < FRAME_SIZE - j; i++){
				first_value = sample_array[i + (80 * k)];
				second_value = sample_array[(i + j) + (80 * k)];
				r[j] = r[j] + ((first_value*hamming_window[i])*(second_value*hamming_window[i + j]));
			}
			//cout << r[j] << endl;
		}
		for (int j = 0; j <= p; j++){
			//printf("R[%d] : %Lf \n", j, r[j]);
			out << fixed << r[j] << " ";
			r[j] = 0;
		}
		out << endl;
		//cout << endl;
	}
	out.close();
	cout << endl;
}

//Calculate Ai's using Levenson-Durbin
void levenson_durbin(long double r[]){
	int i, j;
	long double summation = 0;

	E[0] = r[0];
	for (i = 1; i <= p; i++){
		summation = 0.0;
		for (j = 1; j <= i - 1; j++){
			summation += alpha[j][i - 1] * r[i - j];
		}
		k[i] = (r[i] - summation) / E[i - 1];
		alpha[i][i] = k[i];
		for (j = 1; j <= i - 1; j++){
			alpha[j][i] = alpha[j][i - 1] - (k[i] * alpha[i - j][i - 1]);
		}
		E[i] = (1 - (k[i] * k[i]))*E[i - 1];
	}
	//cout << endl;

	//Calculating the Ai's as the last last column of index i.e. 12th column of matrix alpha
	out.open(ai_file, std::ios_base::app);
	a[0] = 0.0;
	out << a[0] << " ";
	for (int i = 1; i <= p; i++){
		a[i] = alpha[i][12];
		//printf("\nA[%d] : %Lf ", i, a[i]);
		out << fixed << a[i] << " ";
	}
	out << endl;
	out.close();
}

//calling levenson Durbin for 85 frames one after the another
void calculate_Ais(){
	//cout << "\n........Writing Ai's Frames to file.........." << endl;
	int abc = 1;
	string temp;
	int count = 0;
	in.open(ri_file);
	out.open(ai_file);
	out.close();
	for (int k = 0; k < FRAMES; k++){
		while (count <= 12){
			in >> temp;
			r[count] = stod(temp);
			count++;
		}
		levenson_durbin(r);
		count = 0;
	}
	in.close();
}

//Calculate Ci's
void calculate_Cis(){
	//cout << "\n........Writing Ci's for 5 Frames to files.........." << endl;
	string temp;
	long double summation = 0;
	int count = 0;
	long double a1 = 0.0, a2 = 0.0;
	in.open(ai_file);
	in1.open(ri_file);
	out.open(ci_file);
	out.close();

	//Calculating Cis by reading files of Ais and Ris
	out.open(ci_file/*, std::ios_base::app*/);
	for (int k = 0; k < FRAMES; k++){
		while (count <= 12){
			in >> temp;
			//cout << temp << endl;
			a[count] = stod(temp);
			in1 >> temp;
			r[count] = stod(temp);
			count++;
		}
		c[0] = 2 * log(r[0]) / log(2.0);
		for (int m = 1; m <= p; m++){
			summation = 0;
			for (int k = 1; k <= m - 1; k++)
			{
				summation += ((double)k / m)*c[k] * a[m - k];
			}
			c[m] = a[m] + summation;
		}
		for (int m = 1; m <= p; m++){
			out << fixed << c[m] << " ";
			//printf("\nC[%d] : %Lf ", m, c[m]);
			c[m] = 0;
		}
		//cout << endl;
		out << endl;
		count = 0;
	}
	out.close();
	in1.close();
	in.close();
}

//To calculate C Prime(Raised Sine wave)
void calculate_c_prime(){
	//cout << "\n........Writing Ci Prime Frames to file.........." << endl;
	string temp;
	int count = 0;
	in.open(ci_file);
	
	//Calculting CiPrimes by reading Ci values from File
	out1.open(c_prime_file);
	
	for (int k = 0; k < FRAMES; k++){
		while (count < 12){
			in >> temp;
			c[count] = stod(temp);
			count++;
		}
		for (int m = 1; m <= p; m++){
			w[m] = 1.0 + 6.0 * sin((22.0 / 7.0)*m / 12.0);
			c_prime[m] = c[m - 1] * w[m];
		}
		for (int m = 1; m <= p; m++){
			//printf("\nC_Prime[%d] : %Lf ", m, c_prime[m]);
			out << c_prime[m] << " ";
			out1 << c_prime[m] << " ";
			//out1 << c_prime[m] << " ";
			c[m] = 0;
			c_prime[m] = 0;
		}
		//cout << endl;
		out << endl;
		out1 << endl;
		//out1 << endl;
		count = 0;
	}
	in.close();
	out.close();

	//out1.close();
}

void calculate_cepstral_values(){
	

	get_hamming_window();
	calculate_dc_shift();
	calculate_normalization_ratio();
	dc_normalize();

	calculate_Ris();
	calculate_Ais();
	calculate_Cis();
	calculate_c_prime();
}

//TO READ CODEBOOK FROM FILE
void read_code_book(){
	in.open(code_book_file);
	for (int i = 0; i < CODE_BOOK_SIZE; i++){
		for (int j = 0; j < 12; j++){
			in >> code_book[i][j];
		}
	}
	in.close();
}

//TO PRINT CODEBOOK
void print_code_book(){
	cout << "\n**************** Code Book **************" << endl;
	for (int i = 0; i < CODE_BOOK_SIZE; i++){
		for (int j = 0; j < 12; j++){
			cout << " " << code_book[i][j];
		}
		cout << endl;
	}
}

//Calculating Tokhura's Distance Using Code Book
void calculate_tokhura_distance(long double c[12], int index){
	int  min_index = 0;
	long double min = 99999;
	long double sum[CODE_BOOK_SIZE] = { 0 };
	string temp, temp1;

	for (int j = 0; j < CODE_BOOK_SIZE; j++){
		for (int i = 0; i < 12; i++){
			sum[j] += tokhura_weight[i] * (c[i] - code_book[j][i])*(c[i] - code_book[j][i]);
		}
		if (sum[j] < min){
			min = sum[j];
			min_index = j;
		}
	}

	obs[index+1] = min_index + 1;

	//cout << obs[index+1] << " ";
}

void read_ci_values(){
	int i, j;
	string temp;
	in.open(c_prime_file);
	for (i = 0; i < T; i++){
		for (j = 0; j < 12; j++){
			in >> temp;
			c_prime_arr[i][j] = stold(temp);
		}
	}
	in.close();
}

void generate_observation_sequence(){
	read_ci_values();
	int i;
	//cout << "\nObservation Sequence:" << endl;
	for (i = 0; i < FRAMES; i++){
		calculate_tokhura_distance(c_prime_arr[i], i);
	}
	//write_observation_sequence();
}

void write_observation_sequence(){

	int i;
	for (i = 0; i < T; i++)
		out << obs[i] << "\t";
	out << endl;
}

//TO READ A MATRIX FROM FILE
void read_Aij_values(char filename[]){
	in.open(filename);
	string temp;
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			in >> temp;
			a_i_j[i][j] = stold(temp);
		}
	}
	/*cout << "A:" << endl;
	for (int i = 1; i <= N; i++){
	for (int j = 1; j <= N; j++){
	printf("%.30Lf  ", a_i_j[i][j]);
	}
	cout << endl;
	}*/
	in.close();
}

//TO READ B MATRIX FROM FILE
void read_Bij_values(char filename[]){
	in.open(filename);
	string temp;
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= M; j++){
			in >> temp;
			b_i_j[i][j] = stold(temp);
		}
	}
	/*cout << "B:" << endl;
	for (int i = 1; i <= N; i++){
	for (int j = 1; j <= M; j++){
	printf("%.30Lf\n", b_i_j[i][j]);
	}
	cout << endl << endl;
	}*/
	in.close();
}

//TO READ PI VALUES FROM FILE
void read_pi_values(char filename[]){
	in.open(filename);
	string temp;
	int index = 1;
	while (in >> temp){
		pi[index++] = stold(temp);
	}
	in.close();
	/*cout << "Pi values:" << endl;
	for (int i = 1; i < index; i++)
	printf("%.30Lf\n", pi[i]);*/
}


void read_average_model(int iteration){
	char index[3];
	strcpy(a_average_file, "average_models/A_avg_");
	sprintf(index, "%d", iteration);
	strcat(a_average_file, index);
	strcat(a_average_file, ".txt");
	//cout << a_average_file << endl;
	read_Aij_values(a_average_file);

	strcpy(b_average_file, "average_models/B_avg_");
	sprintf(index, "%d", iteration);
	strcat(b_average_file, index);
	strcat(b_average_file, ".txt");
	//cout << b_average_file << endl;
	read_Bij_values(b_average_file);

	strcpy(pi_average_file, "average_models/pi_avg_");
	sprintf(index, "%d", iteration);
	strcat(pi_average_file, index);
	strcat(pi_average_file, ".txt");
	//cout << pi_average_file << endl;
	read_pi_values(pi_average_file);
}

//TO PERFORM THE FORWARD PROCEDURE
void forward_procedure(int iteration){

	int i, j, t;
	long double sum = 0;
	//for (i = 0; i < T; i++)
		//cout << "OBSss:" << obs[i] << " ";

	for (i = 1; i <= N; i++){
		alpha_arr[i][1] = pi[i] * b_i_j[i][obs[1]];		
		//cout << alpha[i][1] << endl;
	}

	//CALCULATING ALPHA MATRIX
	for (t = 1; t <= T - 1; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				
				sum += alpha_arr[i][t] * a_i_j[i][j];
			}

			alpha_arr[j][t + 1] = sum * b_i_j[j][obs[t + 1]];
		}
	}

	//PRINTING ALPHA
	/* << "\nAlpha Matrix:" << endl;
	for (i = 1; i <= N; i++){
	for (j = 1; j <= T; j++){
	printf("%Le\n", alpha[i][j]);
	}
	cout << endl << endl;
	}*/



	//WRITING ALPHA TO FILE
	out.open(alpha_file);
	for (i = 1; i <= T; i++){
		for (j = 1; j <= N; j++){
			out << /*fixed << setprecision(FILE_DECIMAL_PRECISION_VALUE) <<*/ alpha_arr[j][i] << " ";
		}
		out << endl;
	}
	out.close();


	//CALCULATING PROBABILITY OF OBSERVATION SEQUENCE GIVEN MODEL
	sum = 0;
	for (i = 1; i <= N; i++){
		sum += alpha_arr[i][T];
	}

	Pobs_model = sum;
	if (Pobs_model > max_Pobs_model){
		max_Pobs_model = Pobs_model;
		output_model_name = iteration;
	}
	cout << "Digit:"<<iteration<<"\tP(obs/model) : " << Pobs_model <<endl;
}


void solution_to_problem1(int iteration){	
	forward_procedure(iteration);
}

void check_detection(int digit){	
	cout << "\nActual Digit File:" << digit<<endl;
	cout << "Detected Digit:" << output_model_name << endl;
	if (digit == output_model_name){
		count_check++;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	int i, j ,k;
	int utterence_counter = 0;
	int choice = 0;

	//cout << ".....Recording will be for 3 seconds......" << endl;
	//cout << "....... Recording Silence......" << endl;
	//system("Recording_Module.exe 3 silence.wav silence_file.txt");
	//cout << "\nSilence recorded. **Press Enter to record your DIGIT**" << endl;
	//system("Recording_Module.exe 3 input.wav input.txt");
	//cout << "\nRecording successfull **Press ENTER to proceed with program**" << endl;

	cout << "\nSELECT A CHOICE:" << endl;
	cout << "Press 1. To Test using TEST DATA in Test Folder" << endl;
	cout << "Press 2. To Give Your Filename (*File Must be in 'recorded digits' Folder*)" << endl;
	cout << "Enter Your Choice ?" << endl;
	cin >> choice;

	if (choice == 1){
		count_check = 0;
		for (i = 0; i <= DIGITS; i++){
			char index[3];

			for (j = 21; j <= 30; j++){
				utterence_counter++;
				strcpy(input_file, "Test/184101027_");
				sprintf(index, "%d", i);
				strcat(input_file, index);
				strcat(input_file, "_");
				sprintf(index, "%d", j);
				strcat(input_file, index);
				strcat(input_file, ".txt");

				cout << "\n.................................................................................................." << endl;
				cout << "Reading Input from : " << input_file << endl;
				
				calculate_cepstral_values();
				read_code_book();
				//print_code_book();
				generate_observation_sequence();
				out.open(obs_sequence_file);
				write_observation_sequence();
				out.close();
				output_model_name = 0;
				max_Pobs_model = 0;
				for (k = 0; k <= 9; k++){
					read_average_model(k);
					solution_to_problem1(k);
					
				}

				check_detection(i);
				
				
			}
		}
		cout << "\n\nTotal Files Checked:" << utterence_counter << endl;
		cout << "Correct Recognized:" << count_check << endl;
	}
	else if (choice == 2){
		strcpy(input_file, "recorded_digits/");
		cout << "\nEnter File Name:" << endl;
		cin >> custom_input_file;
		strcat(input_file, custom_input_file);
		cout << "\n.................................................................................................." << endl;
		cout << "Reading Input from : " << input_file << endl;
		
		calculate_cepstral_values();
		read_code_book();
		//print_code_book();
		generate_observation_sequence();
		out.open(obs_sequence_file);
		write_observation_sequence();
		out.close();
		output_model_name = 0;
		max_Pobs_model = 0;
		for (k = 0; k <= 9; k++){
			read_average_model(k);
			solution_to_problem1(k);
		}

		cout << "Detected Digit:" << output_model_name << endl;
	}
	else{
		cout << "\nPlease Enter a Valid Choice" << endl;
	}
	return 0;
}

