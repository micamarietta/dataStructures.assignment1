//Mica Marietta
//Data Structures
//Assignment 1
//2318435

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

//calculates and returns the total amount of nucloetides
int calcSum(string totalStrands){
  int sum = 0;
  sum += totalStrands.length();
  return sum;
}

//calculates and returns the mean (avg number of nucleotides per strand)
double calcMean(string totalStrands, double totalLines){
  double mean = (totalStrands.length()/totalLines);
  return mean;
}

//calculates and returns the variance (sum portion of the variance calculated through reading file)
double calcVariance(double varianceSum, double totalLines){
  double variance = (varianceSum/(totalLines-1));
  return variance;
}

//calculates and returns the standard deviation from the variance
 double calcStandDev(double var){
   double standD = sqrt(var);
   return standD;
 }

//calculates and returns the probability of each nucleotide out of all nucleotides
 double calcProb(string allLength, int numofChar){
   double prob = (double) numofChar / (double) allLength.length();
   return prob;
 }

//-----------------------------------------------------------------

int main(int argc, char **argv){
  bool continuePrompting = true;
  bool comLineInput = true;

  while(continuePrompting){
    //declare all variables
    ifstream inFS;
    ofstream outFS;
    string inputFileName = "";
    double lineCount = 0.0;
    string dnaStrand = "";
    string allStrands = "";
    double mean = 0.0;
    int sum = 0;
    double variance = 0.0;
    double standardDeviation = 0.0;

    //probability variables
    int numOfA = 0;
    int numOfT = 0;
    int numOfG = 0;
    int numOfC = 0;
    int numOfAA = 0;
    int numOfTT = 0;
    int numOfGG = 0;
    int numOfCC = 0;
    int numOfAC = 0;
    int numOfAT = 0;
    int numOfAG = 0;
    int numOfTC = 0;
    int numOfTG = 0;
    int numOfGC = 0;

    //gaussian distribution variables
    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;
    string gaussian = "";

    //if this is the user's 1st and on time of inputting take in user input
    //not command line input
    if(!comLineInput){
      cout << "What is the name of the file that you want to read?" << endl;
      cin >> inputFileName;
    }
    else{
      //take in a file from COMMAND LINE input
      if(argc > 1){
        inputFileName = argv[1];
        cout << "Input file name taken in by command line input." << endl;
        comLineInput = false;
      } else{
        cout << "Error: you did not enter the name of the file to read." << endl;
        return 1;
      }
    }

    //open input file
    inFS.open(inputFileName);

  //indicate if the file opens
    if(inFS.is_open()){
      cout << "Input file opened." << endl;
    }

    //-----------------calculating sum and mean--------------------

    //read from the input file
    //concatenate a string of all DNA nucleotides for sum
    //count the amount of strands in the file for mean
    while(!inFS.eof()){
      inFS >> dnaStrand;
      if(!inFS.fail()){
        allStrands += dnaStrand;
        lineCount ++;
      }
    }

    //assign values to global variables
    sum = calcSum(allStrands);
    mean = calcMean(allStrands, lineCount); //need for variance

    //close input file
    inFS.close();

    //-----------------calculating variance--------------------

    //read from file again for variance
    inFS.open(inputFileName);

  //indicate if file was opened
    if(inFS.is_open()){
      cout << "Input file opened again." << endl;
      variance = 0;
    }

    //read each line from file again
    //get length of each line and plug in to variance equation
    while(!inFS.eof()){
      inFS >> dnaStrand;
      if(!inFS.fail()){
        variance += pow((dnaStrand.length() - mean),2); //equation is (lengthofoneline - mean)^2
      }
    }

  //close file again
    inFS.close();


    //------------------calculating probability and bi-gram probabilities------------------------

    //read from file again for probabilities
    inFS.open(inputFileName);

    //indicate if file was opened
    if(inFS.is_open()){
      cout << "Input file opened again." << endl;
    }

    //read each line from file again
    //determine whether each character is an A,T,G, or C and add to probability
    //determine bigram probs (AA,AC,AT,AG,TT,TG,TC,GG,GC,CC)
    while(!inFS.eof()){
      inFS >> dnaStrand;
      if(!inFS.fail()){
        //loop through characters
        for(int i = 0; i < dnaStrand.length();++i){

          if(toupper(dnaStrand[i]) == 'A'){
              numOfA++;

              if(toupper(dnaStrand[i+1]) == 'A'){
                numOfAA++;
              }
              else if(toupper(dnaStrand[i+1]) == 'T'){
                numOfAT++;
              }
              else if(toupper(dnaStrand[i+1]) == 'C'){
                numOfAC++;
              }
              else if(toupper(dnaStrand[i+1]) == 'G'){
                numOfAG++;
              }
          }
          else if(toupper(dnaStrand[i]) == 'T'){
              numOfT++;

              if (toupper(dnaStrand[i+1]) == 'T'){
                numOfTT++;
              }
              else if(toupper(dnaStrand[i+1]) == 'C'){
                numOfTC++;
              }
              else if(toupper(dnaStrand[i+1]) == 'G'){
                numOfTG++;
              }
          }
          else if(toupper(dnaStrand[i]) == 'G'){
              numOfG++;

              if(toupper(dnaStrand[i+1]) == 'C'){
                numOfGG++;
              }
              else if(toupper(dnaStrand[i+1]) == 'G'){
                numOfGC++;
              }
          }
          else if(toeupper(dnaStrand[i]) == 'C'){
              numOfC++;

              if(toupper(dnaStrand[i+1]) == 'C'){
                numOfCC++;
              }
          }
        }
      }
    }

    //open output file
    //append to file
    outFS.open("micamarietta.out", ios::app);

    //indicate if output file was opened
    if(outFS.is_open()){
      cout << "micamarietta.out opened." << endl;
    }

    //assign global variables to actual values through calculations
    sum = calcSum(allStrands);
    mean = calcMean(allStrands, lineCount);
    variance = calcVariance(variance, lineCount);
    standardDeviation = calcStandDev(variance);

    //write statements to out file
    outFS << "Mica Marietta" << endl;
    outFS << "Data Structures" << endl;
    outFS << "2318435" << endl;
    outFS << "Sum:   " << sum << endl;
    outFS << "Mean:   " << mean << endl;
    outFS << "Variance:   " << variance << endl;
    outFS << "Standard Deviation:   " << standardDeviation << endl;
    outFS << "Probability of A:    " << calcProb(allStrands, numOfA) << endl;
    outFS << "Probability of T:    " << calcProb(allStrands, numOfT) << endl;
    outFS << "Probability of G:    " << calcProb(allStrands, numOfG) << endl;
    outFS << "Probability of C:    " << calcProb(allStrands, numOfC) << endl;
    outFS << "Probability of AA:   " << calcProb(allStrands, numOfAA) << endl;
    outFS << "Probability of AT:   " << calcProb(allStrands, numOfAT) << endl;
    outFS << "Probability of AG:   " << calcProb(allStrands, numOfAG) << endl;
    outFS << "Probability of AC:   " << calcProb(allStrands, numOfAC) << endl;
    outFS << "Probability of TT:   " << calcProb(allStrands, numOfTT) << endl;
    outFS << "Probability of TG:   " << calcProb(allStrands, numOfTG) << endl;
    outFS << "Probability of TC:   " << calcProb(allStrands, numOfTC) << endl;
    outFS << "Probability of GG:   " << calcProb(allStrands, numOfGG) << endl;
    outFS << "Probability of GC:   " << calcProb(allStrands, numOfGC) << endl;
    outFS << "Probability of CC:   " << calcProb(allStrands, numOfCC) << endl;
    outFS << endl;

    //produce 1000 lines that follow the gaussian distribution
    for(int i = 0; i < 1000; ++i){

      //generate random floats for a and b
      a = (rand()/ (double)(RAND_MAX));
      b = (rand()/ (double)(RAND_MAX));

      //calculate standard gaussian and length of new strings
      c = abs((sqrt(-2 * log(a))) * cos (2 * M_PI * b));
      d = (double)standardDeviation * c + mean;
      d = round(d);

      //calculate probabilites for outputting 1000 strings
      double tempA = d * calcProb(allStrands, numOfA);
      double tempT = d * calcProb(allStrands, numOfT);
      double tempC = d * calcProb(allStrands, numOfC);
      double tempG = d * calcProb(allStrands, numOfG);

      //generate each nucleotide based on gaussian probability
      for(int j = 0; j < tempA; ++j){
        gaussian += 'A';
      }
      for(int j = 0; j < tempT; ++j){
        gaussian += 'T';
      }
      for(int j = 0; j < tempC; ++j){
        gaussian += 'C';
      }
      for(int j = 0; j < tempG; ++j){
        gaussian += 'G';
      }

      //output string to .out file and reset it
      outFS << gaussian << endl;
      gaussian = "";
    }

    outFS.close();

    //ask the user if they would like to continue
    cout << "Would you like to continue? (answer yes or no)" << endl;
    string answer = "";
    cin >> answer;
    if(answer == "yes"){
      continuePrompting = true;
    } else{
      continuePrompting = false;
      cout << "Ending program." << endl;
    }
  }

  return 0;
}
