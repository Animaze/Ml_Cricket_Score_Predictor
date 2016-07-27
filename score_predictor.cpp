#include <bits/stdc++.h> 
using namespace std; 
#include<stdlib.h>
//player_info gives the batting and bowling stats of the player.
#define MAX 700
#define INT_MAX 1000000009
#define N 4
int total_match=100;
double input[MAX][MAX],output[MAX][MAX],inputTranspose[MAX][MAX],productA[N][N],f,y[N][N],score_matrix[MAX][2],wicket_matrix[MAX][2],productB[MAX][MAX],weight_matrix[MAX][2],final[MAX][MAX],identity[MAX][MAX];   //product A is the multiplication of input & inputTranspose
int inputRow=2*total_match,inputCol =N,outputRow=N,outputCol=1;
double lamda=0,error,min_error=INT_MAX,min_lamda=0;                      
double s_inv[N][N];
double productLA[MAX][MAX],actualScores[MAX],calcScores[MAX],validate[MAX][N],score_val[MAX],min_outSampleErr=INT_MAX,outSampleErr[5000];     //outSampleError matrix 1 error for a particular lamda.
double D[N][N],C[N][N];
int set_assign[1000][1000],total_values_in_set[1000],set_number,set_of[1000];
double mean_eout=0.0;
int lamda_ctr,score_total=0;
int validation_folds=10;
class player_info{
		 
	 public:
     int player_matches,odi_no; 
     char name[40],country[20];
     int bat_inns,bat_NO,bat_runs,bat_HS,bat_BF,bat_100,bat_50;
     double bat_SR,bat_avg; 
	 int bowl_inns,bowl_wkt,bowl_runs,bowl_balls,bowl_5W,bowl_4W;              //player_cap and the country would be the primary key.
	 double bowl_SR,bowl_avg,bowl_economy,bowl_overs;
	 double sum_avg,over_per_match;	                                                        //sum of bat and bowl avgs. higher->BATSMEN
};
class match{
	public:
	int odi_no,d_l,skew;
	char teamA_name[30],teamB_name[30];
	char ground[20],H_N,winner,toss,bat_first;
	double overA,overB,maxover_A,maxover_B;
	int wktA,wktB,scoreA,scoreB;
	double teamA_bat_avg,teamB_bat_avg,teamA_bowl_avg,teamB_bowl_avg;                    //These gives the batting and bowling abilities of the two teams.
	double teamA_bat_SR,teamB_bat_SR,teamA_bowl_SR,teamB_bowl_SR;
	double teamA_bowl_eco,teamB_bowl_eco;
};
class team_inp{
	public:
	int d_l;
	int H_N,winner,toss,bat_first;
	double maxover;
	double team_bat_avg,opp_bowl_avg;                    //These gives the batting and bowling abilities of the two teams.
	double team_bat_SR,opp_bowl_SR;
	double opp_bowl_eco;
};

class score_out{
	public:
		int score;
};
class wicket_out{
	public:
		int wicket;
};


//Linear Decomposition of (x'*x)
void LU(double(*D)[N][N],int n){
	int i,j,k,m;
	double x;
    for(k=0;k<=n-1;k++){
		for(j=k+1;j<=n;j++){
		    x=(*D)[j][k]/(*D)[k][k];
		    for(i=k;i<=n;i++){  
		       (*D)[j][i]=(*D)[j][i]-x*(*D)[k][i];
	        }
		    (*D)[j][k]=x;
		}
    }
}

//Calculates the inverse of (x'*x + lamda*I)
int calc_inverse(){
	int i,j,n,m;
    double d[N];
    double x,y[N];
    n=N-1;

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		       D[i][j]=productLA[i][j];
		}
	}	
	
    for(m=0;m<=N;m++){ 
		for(j=0;j<=N;j++)
	  		C[m][j]=D[m][j];
    }
       
    LU(&D,n);

    for(m=0;m<N;m++){ 
    
    	for(i=0;i<N;i++)
        	d[i]=0.0;
		d[m]=1.0;
		for(i=0;i<=n;i++){
			 x=0.0; 
		  	for(j=0;j<=i-1;j++){
			   x=x+D[i][j]*y[j];
			}
	 	  	y[i]=(d[i]-x);
		}

		for(i=n;i>=0;i--){ 
			x=0.0; 
		  	for(j=i+1;j<=n;j++){
			  	 x=x+D[i][j]*s_inv[j][m];
			}
	 		s_inv[i][m]=(y[i]-x)/D[i][i];
		}
	}
	FILE *ptr_matrix;
	ptr_matrix =fopen("Matrix_after_inversion.txt", "w");

	 for(m=0;m<N;m++){ 
        for(i=0;i<N;i++)
        	fprintf(ptr_matrix,"%lf ",s_inv[m][i]);
	}
	fclose(ptr_matrix);

/*	Checks the correctness of inverse.
printf("The product\n");
        for(m=0;m<N;m++)
        { 
	 for(j=0;j<N;j++)
	  {
           x=0.0; 
	   for(i=0;i<N;i++)
	    { x=x+C[m][i]*s_inv[i][j]; }
          printf("%lf ",x ); 
	}
	printf("\n");
	
}*/
	return 0;
}

int compare_bowl(const void *a,const void *b){
	return((*(player_info*)b).over_per_match-(*(player_info*)a).over_per_match);
}
int compare(const void *a,const void *b){
	return((*(player_info*)b).bat_avg-(*(player_info*)a).bat_avg);
}

class player_info p[300];
class match m[350];
class team_inp t[700];
class score_out s[700];
class wicket_out w[700];



void readValuesIntoMatrix(int val){                          //Reads values into the input matrix and also forms the input transpose matrix
	int matchCtr = 0,j1 = 0, j;
	if(val==0){
		matchCtr = 0;
		inputRow=2*total_match;
		for(int i=0;i<inputRow;i++){
				input[i][0]=1;
				input[i][1]=t[matchCtr].team_bat_avg;
				input[i][2]=t[matchCtr].opp_bowl_avg;
				input[i][3]=t[matchCtr].team_bat_SR;
				//input[i][4]=t[matchCtr].opp_bowl_SR;
				//input[i][4]=t[matchCtr].H_N;
			//	input[i][4]=t[matchCtr].opp_bowl_eco;			
		    	matchCtr++;	
		
		}
		for(int i=0;i<inputRow;i++){
			score_matrix[i][0] = s[i].score;
		}
		
		for(int i=0;i<inputRow;i++){
			for(int j=0;j<inputCol;j++){
				inputTranspose[j][i] = input[i][j];
			}
		}
	}
	else{   
	        j1=0;
			j=0;
			inputRow=2*total_match;
	        matchCtr = 0;
			for(int i=0;i<inputRow;i++){
			    if(set_of[i]==set_number){
 				validate[j1][0]=1;
				validate[j1][1]=t[matchCtr].team_bat_avg;
				validate[j1][2]=t[matchCtr].opp_bowl_avg;
				validate[j1][3]=t[matchCtr].team_bat_SR;
				validate[j1][4]=t[matchCtr].opp_bowl_eco;
				j1++;			
		    	matchCtr++;	               
			    continue;}
				input[j][0]=1;
				input[j][1]=t[matchCtr].team_bat_avg;
				input[j][2]=t[matchCtr].opp_bowl_avg;
				input[j][3]=t[matchCtr].team_bat_SR;
				input[j][4]=t[matchCtr].opp_bowl_eco;			
		    	matchCtr++;	
		    	j++;
		   }
		    
		    j=0;
		    j1=0;
			for(int i=0;i<inputRow;i++){
				if(set_of[i]==set_number){
					score_val[j1]=s[i].score;
					j1++;
				    continue;
					}
				score_matrix[j][0] = s[i].score;
				j++;
			}
			inputRow=j;
			for(int i=0;i<inputRow;i++){
				for(int k=0;k<inputCol;k++){
					inputTranspose[k][i] = input[i][k];
				}
			}
		  	

	
		
	}

	
}

//ProductA gives the product of x' & x
void firstMultiply(){  

		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputCol;j++){
				productA[i][j] = 0;
			}
		}	
         //Multiply input transpose and input matrix , x'*x
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputCol;j++){
				for(int k=0;k<inputRow;k++){
					productA[i][j] += inputTranspose[i][k]*input[k][j]; 
				}
			}
		}
}

//Adds Identity(lamda *I) Matrix to ProductA
void firstAdd(){
	
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputCol;j++){
				productLA[i][j] = 0;
			}
		}	
	
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputCol;j++){
				if(i==j)
					identity[i][j]=1;
				else
					identity[i][j]=0;
			productLA[i][j] = productA[i][j] + lamda*identity[i][j];}
		}		
}

//Multiplies the inverse with the input TRANSPOSE matrix, i.e , ((x'*x)^-1)*x'
void secondMultiply(){											
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputRow;j++){
				productB[i][j] = 0;
			}
		}
		
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<inputRow;j++){	 
				for(int k=0;k<inputCol;k++){
					productB[i][j] += y[i][k]*inputTranspose[k][j]; 
				}
			}
		}
		
return;
}

//Multiplies the pseudo-inverse with the output matrix, i.e , (((x'*x)^-1)*x')  *  y
void thirdMultiply(){	
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<outputCol;j++){
				weight_matrix[i][j] = 0;
			}
		}
										
		for(int i=0;i<inputCol;i++){
			for(int j=0;j<outputCol;j++){
				for(int k=0;k<inputRow;k++){
					weight_matrix[i][j] += productB[i][k]*score_matrix[k][j]; 
				}
			}
		}
return;
}

//Calculates the predicted scores	
void calcResult(){
	
		for(int i=0;i<inputRow;i++){
			for(int j=0;j<outputCol;j++){
				final[i][j]=0;
			}
		}	
		
		for(int i=0;i<inputRow;i++){
			for(int j=0;j<outputCol;j++){
				for(int k=0;k<inputCol;k++){
					final[i][j] += input[i][k]*weight_matrix[k][j];
				}
			}
		}	
}	

//Calcutes Squared Error
int calcError(){

	for(int i=0;i<inputRow;i++){	
		actualScores[i]=s[i].score;
	}
	
	error = 0;
	for(int i=0;i<inputRow;i++){
		error += (actualScores[i] - final[i][0])*(actualScores[i] - final[i][0]); 
	}
	
	if(error < min_error ){
		min_error = error;
		min_lamda = lamda;
	}
}

//Divides the inputs into training and validation set
void divide_training_set(){
	int i,set,total_sets,x;
	total_sets=validation_folds;
	for(i=0;i<total_sets;i++)
	   total_values_in_set[i]=0;
	   
	for(i=0;i<2*total_match;i++){
	      set=(rand()+9)%total_sets;
	      while(total_values_in_set[set]>=(2*total_match)/validation_folds+1){
	      	set=(rand()+9)%total_sets;
		  }
		  set_assign[set][total_values_in_set[set]++]=i;
		  set_of[i]=set;
	}
}
void calcResultForValidation(){
	
		for(int i=0;i<total_values_in_set[set_number];i++){
			for(int j=0;j<outputCol;j++){
				final[i][j]=0;
			}
		}	
		
		for(int i=0;i<total_values_in_set[set_number];i++){
			for(int j=0;j<outputCol;j++){
				for(int k=0;k<inputCol;k++){
					final[i][j] += validate[i][k]*weight_matrix[k][j];
				}
			}
		}
		
}	
void calcErrorForValidation(){
	
	calcResultForValidation();
	
	for(int i=0;i<total_values_in_set[set_number];i++){
		outSampleErr[lamda_ctr] += (final[i][0]-score_val[i])*(final[i][0]-score_val[i]);

	}
	


	
}

int main()
	{
		FILE *ptr_player,*ptr_match,*ptr_score,*ptr_inp;

		ptr_player =fopen("player.txt", "r");
		ptr_match =fopen("match.txt", "r");
	    
	    ptr_score = fopen("sc1.txt","w");
		char tempc[50];
		if (!ptr_player)
			return 1;
		if(!ptr_match)
		 return 1;
		if(!ptr_score)
		 return 1;
			
		int i=0,j,match_ctr=total_match,x,temp,i1;
	
		while(match_ctr--){
			
			fscanf(ptr_match,"%d %d %lf %d %d %lf %d %s %c %c %d %c %c %d %lf %lf",&m[i].odi_no,&m[i].wktA,&m[i].overA,&m[i].scoreA,&m[i].wktB,&m[i].overB,&m[i].scoreB,m[i].ground,&m[i].H_N,&m[i].winner,&m[i].d_l,&m[i].toss,&m[i].bat_first,&m[i].skew,&m[i].maxover_A,&m[i].maxover_B);
			x=11;
			j=2;
		    	if(m[i].skew!=0||m[i].d_l!=0||m[i].maxover_A!=50||m[i].maxover_B!=50||m[i].scoreA<160||m[i].scoreB<160){
		    		fscanf(ptr_player,"%s",m[i].teamA_name);
		    		for(i1=0;i1<11;i1++)
	               	fscanf(ptr_player,"%s %d %d %d %d %d %lf %lf %d %d %d %d %d",p[i1].name,&p[i1].player_matches,&p[i1].bat_inns,&p[i1].bat_NO,&p[i1].bat_runs,&p[i1].bat_HS,&p[i1].bat_avg,&p[i1].bat_SR,&p[i1].bat_100,&p[i1].bat_50,&temp,&temp,&temp);
			        fscanf(ptr_player,"%s",m[i].teamA_name);
		    		for(i1=0;i1<11;i1++)
					fscanf(ptr_player,"%s %d %d %lf %d %d %d %lf %lf %lf %d %d %d %d",tempc,&temp,&p[i1].bowl_inns,&p[i1].bowl_overs,&temp,&p[i1].bowl_runs,&p[i1].bowl_wkt,&p[i1].bowl_avg,&p[i1].bowl_economy,&p[i1].bowl_SR,&temp,&temp,&temp,&temp);
					fscanf(ptr_player,"%s",m[i].teamB_name);
		    		for(i1=0;i1<11;i1++)
					fscanf(ptr_player,"%s %d %d %d %d %d %lf %lf %d %d %d %d %d",p[i1].name,&p[i1].player_matches,&p[i1].bat_inns,&p[i1].bat_NO,&p[i1].bat_runs,&p[i1].bat_HS,&p[i1].bat_avg,&p[i1].bat_SR,&p[i1].bat_100,&p[i1].bat_50,&temp,&temp,&temp);
			        fscanf(ptr_player,"%s",m[i].teamB_name);
		    		for(i1=0;i1<11;i1++) 
					fscanf(ptr_player,"%s %d %d %lf %d %d %d %lf %lf %lf %d %d %d %d",tempc,&temp,&p[i1].bowl_inns,&p[i1].bowl_overs,&temp,&p[i1].bowl_runs,&p[i1].bowl_wkt,&p[i1].bowl_avg,&p[i1].bowl_economy,&p[i1].bowl_SR,&temp,&temp,&temp,&temp);
				  	continue;
				 }

			fprintf(ptr_score,"%d\n%d\n",m[i].scoreA,m[i].scoreB);
			while(j--){
				
				i1=0;
				if(j==1)
					fscanf(ptr_player,"%s",m[i].teamA_name);
				else
					fscanf(ptr_player,"%s",m[i].teamB_name);
				x=11;
				while(x--){
					p[i1].sum_avg=0;
					fscanf(ptr_player,"%s %d %d %d %d %d %lf %lf %d %d %d %d %d",p[i1].name,&p[i1].player_matches,&p[i1].bat_inns,&p[i1].bat_NO,&p[i1].bat_runs,&p[i1].bat_HS,&p[i1].bat_avg,&p[i1].bat_SR,&p[i1].bat_100,&p[i1].bat_50,&temp,&temp,&temp);	
					p[i1].sum_avg=p[i1].bat_avg;
					i1++;
					}
		
				if(j==1)
					fscanf(ptr_player,"%s",m[i].teamA_name);
				else
					fscanf(ptr_player,"%s",m[i].teamB_name);
				x=11;
				i1=0;
				while(x--){
					fscanf(ptr_player,"%s %d %d %lf %d %d %d %lf %lf %lf %d %d %d %d",tempc,&temp,&p[i1].bowl_inns,&p[i1].bowl_overs,&temp,&p[i1].bowl_runs,&p[i1].bowl_wkt,&p[i1].bowl_avg,&p[i1].bowl_economy,&p[i1].bowl_SR,&temp,&temp,&temp,&temp);
						
						if(p[i1].bowl_avg==0)
						p[i1].sum_avg=p[i1].sum_avg+100.0;
					//	p[i1].sum_avg=p[i1].sum_avg+p[i1].bowl_avg;
						if(p[i1].player_matches!=0)
						p[i1].over_per_match=p[i1].bowl_overs/(double)(p[i1].player_matches);
						else
						p[i1].over_per_match=0;
						
						i1++;
					}
				if(j==1){
					m[i].teamA_bowl_eco=0; m[i].teamA_bat_avg=0; m[i].teamA_bat_SR=0,m[i].teamA_bowl_avg=0,m[i].teamA_bowl_SR=0;}
				else{
				   	m[i].teamB_bowl_eco=0;m[i].teamB_bat_avg=0;m[i].teamB_bat_SR=0,m[i].teamB_bowl_avg=0,m[i].teamB_bowl_SR=0;}
				    	
				//Sort according to batting measure. In decreasing order of batting ability.
				qsort(p,11,sizeof(player_info),compare);
				double total_exp_bats=0,total_exp_bowl_overs=0;
				for(i1=0;i1<7;i1++){
					if(p[i1].bat_inns>5&&p[i1].bat_avg<70)
					{
						if(j==1){
							 m[i].teamA_bat_avg=m[i].teamA_bat_avg+p[i1].bat_avg;
							 m[i].teamA_bat_SR=m[i].teamA_bat_SR+p[i1].bat_SR;	
							}
						else{
						 	m[i].teamB_bat_avg=m[i].teamB_bat_avg+p[i1].bat_avg;
						 	m[i].teamB_bat_SR=m[i].teamB_bat_SR+p[i1].bat_SR;
						}
						total_exp_bats++;
					}
				}
				
				if(total_exp_bats==0){
					if(j==1){
					   m[i].teamA_bat_avg=30.0;
					   m[i].teamA_bat_SR=75.0;
					}
					else{
						m[i].teamB_bat_avg=30.0;
					    m[i].teamB_bat_SR=75.0;
					}
				}
				else{
					if(j==1){
					   m[i].teamA_bat_avg=m[i].teamA_bat_avg/total_exp_bats;
					   m[i].teamA_bat_SR=m[i].teamA_bat_SR/total_exp_bats;
					}
					else{
						m[i].teamB_bat_avg=m[i].teamB_bat_avg/total_exp_bats;
					    m[i].teamB_bat_SR=m[i].teamB_bat_SR/total_exp_bats;
					}
			    }
			    
				//Sort according to bowling measure. In decreasing order of bowling ability.
			    qsort(p,11,sizeof(player_info),compare_bowl);
			    for(i1=0;i1<6;i1++){
			    	if(p[i1].player_matches>8&&p[i1].over_per_match>1.0)
			    	{
			    		if(j==1){
						 m[i].teamA_bowl_avg=m[i].teamA_bowl_avg+(p[i1].bowl_avg)*(p[i1].over_per_match);
						 m[i].teamA_bowl_SR=m[i].teamA_bowl_SR+p[i1].bowl_SR*(p[i1].over_per_match);
						  m[i].teamA_bowl_eco=m[i].teamA_bowl_eco+p[i1].bowl_economy*(p[i1].over_per_match);	
						}
						else{
						 m[i].teamB_bowl_avg=m[i].teamB_bowl_avg+(p[i1].bowl_avg)*(p[i1].over_per_match);
						 m[i].teamB_bowl_SR=m[i].teamB_bowl_SR+p[i1].bowl_SR*(p[i1].over_per_match);
						  m[i].teamB_bowl_eco=m[i].teamB_bowl_eco+p[i1].bowl_economy*(p[i1].over_per_match);	
						}
			    		total_exp_bowl_overs=total_exp_bowl_overs+p[i1].over_per_match;
					}
				}
				
				if(total_exp_bowl_overs==0){
					if(j==1){
					   m[i].teamA_bowl_avg=45.0;
					   m[i].teamA_bowl_SR=45.0;
					   m[i].teamA_bowl_eco=6.0;
					}
					else{
					   m[i].teamB_bowl_avg=45.0;
					   m[i].teamB_bowl_SR=45.0;
					   m[i].teamB_bowl_eco=6.0;
					}
				}
				else{
					if(j==1){
					   m[i].teamA_bowl_avg=m[i].teamA_bowl_avg/total_exp_bowl_overs;
					   m[i].teamA_bowl_SR=m[i].teamA_bowl_SR/total_exp_bowl_overs;
					   m[i].teamA_bowl_eco=m[i].teamA_bowl_eco/total_exp_bowl_overs;
					}
					else{
					   m[i].teamB_bowl_avg=m[i].teamB_bowl_avg/total_exp_bowl_overs;
					   m[i].teamB_bowl_SR=m[i].teamB_bowl_SR/total_exp_bowl_overs;
					   m[i].teamB_bowl_eco=m[i].teamB_bowl_eco/total_exp_bowl_overs;
					}
			    }
			
			}		
			i++;	
		}
        
        
		total_match=i;
		inputRow=2*total_match;
	    
		divide_training_set();
		fclose(ptr_player);
		fclose(ptr_match);
		fclose(ptr_score);
		
		ptr_inp = fopen("inp.txt","w");
        
		// Construction of input and output matrix.
		   int it=0;
		   for(i=0;i<total_match;i++){
				t[it].d_l=m[i].d_l;
				
		   		if(m[i].H_N=='A')
				   t[it].H_N=1;
				else
				 	t[it].H_N=0;
				 	
				if(m[i].winner=='A')
					t[it].winner = 1;
				else if(m[i].winner=='B')
					t[it].winner = -1;
				else 
					t[it].winner = 0;
				
		   	 	if(m[i].toss == 'A')
					t[it].toss = 1;
				else if(m[i].toss == 'B')
					t[it].toss = -1;
						
				if(m[i].bat_first == 'A')
					t[it].bat_first = 1;
				else if(m[i].bat_first == 'B')
					t[it].bat_first = -1;
					
				 t[it].maxover=m[i].maxover_A;
			   	 t[it].team_bat_avg=m[i].teamA_bat_avg;
			   	 t[it].opp_bowl_avg=m[i].teamB_bowl_avg;
			   	 t[it].team_bat_SR=m[i].teamA_bat_SR;
			   	 t[it].opp_bowl_SR=m[i].teamB_bowl_SR;
			   	 t[it].opp_bowl_eco=m[i].teamB_bowl_eco;
		   	     s[it].score=m[i].scoreA;
		   	     w[it].wicket=m[i].wktA;
		   	     
		   	     fprintf(ptr_inp,"1.0 %lf %lf %lf\n",t[it].team_bat_avg,t[it].opp_bowl_avg,t[it].team_bat_SR);
		   	     
		   	    it++;
		   	  
		   	  	t[it].d_l=m[i].d_l;
		   	  	
		   		if(m[i].H_N=='A')
				   t[it].H_N=-1;
				else 
					t[it].H_N=0;
					
				if(m[i].winner=='A')
					t[it].winner = -1;
				else if(m[i].winner=='B')
					t[it].winner = 1;
				else 
				    t[it].winner = 0;
					
		   	 	if(m[i].toss == 'A')
					t[it].toss = -1;
				else if(m[i].toss == 'B')
					t[it].toss = 1;
					
				if(m[i].bat_first == 'A')
					t[it].bat_first = -1;
				else if(m[i].bat_first == 'B')
					t[it].bat_first = 1;
					
				 t[it].maxover=m[i].maxover_B;
			   	 t[it].team_bat_avg=m[i].teamB_bat_avg;
			   	 t[it].opp_bowl_avg=m[i].teamA_bowl_avg;
			   	 t[it].team_bat_SR=m[i].teamB_bat_SR;
			   	 t[it].opp_bowl_SR=m[i].teamA_bowl_SR;
			   	 t[it].opp_bowl_eco=m[i].teamA_bowl_eco;
			   	 s[it].score=m[i].scoreB;
		   	     w[it].wicket=m[i].wktB;
		   	   	 fprintf(ptr_inp,"1.0 %lf %lf %lf\n",t[it].team_bat_avg,t[it].opp_bowl_avg,t[it].team_bat_SR);
		   	   	 
		   	 	 it++;
		   }
		   
		    fclose(ptr_inp);
		     
		    												
			/*FILE *ptr_matrix;
		    ptr_matrix =fopen("Matrix_to_invert.txt", "w");
	        if (!ptr_matrix)
			return 1;*/
			
			//Multiply x' and x
																
		
	         /*for(int i=0;i<inputCol;i++){
					for(int j=0;j<inputCol;j++){
					fprintf(ptr_matrix,"%lf ",productA[i][j]);
					}
					cout<<endl;	
				}
			fclose(ptr_matrix);*/
		   
		   //Regularization
		   
		 for(set_number=0;set_number<validation_folds;set_number++){

		 	readValuesIntoMatrix(1);
		    firstMultiply();
			for(int i1=0;i1<100;i1++){
				lamda_ctr=i1;
				lamda = (3.985+0.0001*i1); 
				firstAdd();
				calc_inverse();
				for(i=0;i<N;i++){
					for(j=0;j<N;j++){
						y[i][j]=s_inv[i][j];
					}
				}
	    		secondMultiply();
				thirdMultiply();		
				calcErrorForValidation();
				
			}

		}	
	
		
		
		
		  	inputRow=2*total_match;
		  		for(int i1=0;i1<100;i1++){
	            outSampleErr[i1]= outSampleErr[i1]/(double)(2*total_match);
			   outSampleErr[i1]=pow(outSampleErr[i1],0.5);
			   
			   if(outSampleErr[i1]<min_outSampleErr){
			   min_outSampleErr=outSampleErr[i1];
			   min_lamda=3.985+0.0001*i1;
			  // cout<<"hi";
			   }
		  // cout<<outSampleErr[i1]<<endl;
		}
		
		
		    readValuesIntoMatrix(0);
		   /*  cout<<"INPUT"<<endl;
		    for(i=0;i<3;i++){
		    	for(j=0;j<N;j++)
		    	    cout<<inputTranspose[i][j]<<" ";
		    	    cout<<endl;
			}*/
		    firstMultiply();
			lamda = min_lamda;
			firstAdd();
			calc_inverse();
			for(i=0;i<2*total_match;i++)
			score_total+=s[i].score;
			double E_OUT_PERC;
	    	E_OUT_PERC=min_outSampleErr/(double)score_total;
	    	E_OUT_PERC=E_OUT_PERC*(2*total_match);
			E_OUT_PERC=E_OUT_PERC*100.00;
			cout<<"Out of sample err %  : "<<E_OUT_PERC<<"% ";//<<min_outSampleErr<<endl;
	     	//Y matrix holds the inverse of (x'*x)
 			for(i=0;i<N;i++){
				for(j=0;j<N;j++){
					y[i][j]=s_inv[i][j];
				}
			}
    		secondMultiply();
			thirdMultiply();
			calcResult();
			
			cout<<"WEIGHT MATRIX"<<endl;
			for(int i=0;i<outputRow;i++){
				for(int j=0;j<outputCol;j++){
					cout<<weight_matrix[i][j]<<" ";
				}
			cout<<endl;	
			}
			  
            ptr_inp =fopen("inp.txt", "r");  
			for(i=0;i<2*total_match;i++){
				for(j=0;j<N;j++)
					fscanf(ptr_inp,"%lf",&input[i][j]);	
			} 
			fclose(ptr_inp);
			
    	
		
		char check_win;
		int ans=0;

		 cout<<"ANSWER MATRIX\n";
		 double e=0,end=0;
				for(int i=0;i<2*total_match;i++){
					for(int j=0;j<1;j++){
						cout<<final[i][j];
						
						if(s[i].score > final[i][j])
							end = end + s[i].score - final[i][j];
						else
							end = end + final[i][j] - s[i].score;
						
						if(i%2==1){
							if(final[i][j]<=final[i-1][j])
								check_win='A';
							else
								check_win='B';
	
						if(check_win==m[i/2].winner)
							ans++;
						}
					}	
					cout<<endl;	
				}
				cout<<"No of matches in which our prediction(who wins,out of 75 matches) is right :"<<ans<<"\n"<<"In Sample Error % :"<<((end)/(double)score_total)*100.00<<"%";
		
		return  0;
	}






