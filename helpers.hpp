
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <unordered_map>

#ifndef HELPERS
#define HELPERS

// data allocation, contiguous
float** allocate(unsigned int height, unsigned int width, const float& val = 0) {    
    float** ptr = new float*[height]; 
    float* mem = new float[height*width]{ val }; 

    for (unsigned int i = 0; i < height; ++i, mem += width)
        ptr[i] = mem;

    return ptr;
}

void deallocate(float** data) {
   delete [] data[0];  
   delete [] data;     
}

struct SequenceInfo {
    std::vector<char> X, Y; // input sequences
    float match_score = 1.0, mismatch_score = -1.0, gap_penalty = -2.0; // default scoring scheme
    std::vector<char> X_aligned, Y_aligned; // aligned sequences
    
    int rows=0, cols=0, SUB_size=0;// helpers
    int similarity_score = 0, identity_score = 0, gap_count = 0; // output statistics

    // interfaces
    unsigned long gpsa_sequential(float** S);
    unsigned long gpsa_taskloop(float** S, long grain_size, int block_size_x, int block_size_y);
    unsigned long gpsa_tasks(float** S, long grain_size, int block_size_x, int block_size_y);

    SequenceInfo(std::string X_filename, std::string Y_filename) {
        X = load_sequence(X_filename);
        Y = load_sequence(Y_filename);
        rows = X.size()+1;
        cols = Y.size()+1;

        scoring_scheme(1.0, -1.0, -2.0);
    }

    // Traceback, and write aligned sequences
    void traceback_and_save(std::string filename, float** S, bool print=false) {
        std::remove(filename.c_str());

        int i = X.size();
        int j = Y.size();

        // there are multiple solution with a similar score (the output depends on the order of checks)
        while (i > 0 || j > 0) {
            if (i > 0  && S[i][j] == (S[i - 1][j] + gap_penalty)) {
                // left
                X_aligned.insert(X_aligned.begin(),  X[i - 1]);
                Y_aligned.insert(Y_aligned.begin(),  '-');
                gap_count++;
                
                i--;
            } else if (j > 0  && S[i][j] == (S[i][j-1] + gap_penalty)) {
                // up
                X_aligned.insert(X_aligned.begin(),  '-');
                Y_aligned.insert(Y_aligned.begin(),  Y[j - 1]);
                gap_count++;
                
                j--;
            } else {
                //if (i > 0 && j > 0  && (S[i][j] == S[i - 1][j - 1] + (X[i - 1] == Y[j - 1] ? match_score : mismatch_score))) {
                // diagonal top-left
                X_aligned.insert(X_aligned.begin(),  X[i - 1]);
                Y_aligned.insert(Y_aligned.begin(),  Y[j - 1]);
                
                if ((X[i - 1] == Y[j - 1] ? match_score : mismatch_score) == match_score) {
                    similarity_score += 1;
                    if (X[i - 1] == Y[j - 1])
                        identity_score += 1;
                }
                
                i--; j--;
            }
            
        }

        if (print) {
            for ( auto& el: X_aligned)
                std::cout << el;
            std::cout << std::endl;
            for ( auto& el: Y_aligned)
                std::cout << el;
            std::cout << std::endl;
        }

        // save to disk
        std::ofstream ofs(filename, std::ofstream::trunc);
        for ( auto& el: X_aligned) ofs << el;
            ofs << std::endl;
        for ( auto& el: Y_aligned) ofs << el;
            ofs << std::endl;
        ofs.close();
    }

    // Load sequences from input files 
    std::vector<char> load_sequence(std::string filename) {
        std::ifstream ifs(filename);

        std::vector<char> res;

        if (!ifs.good()) {
            std::cerr << "[error]: could not open input file '" << filename << "'!" << std::endl;
            exit(-1);
        }

        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        std::string line;
        while (std::getline(ifs, line)) {
            // str += line; 
            for ( auto& c: line ) {
                res.push_back(c);
            }
        }

        ifs.close();

        return std::move(res);
    }

    // Reset between the runs
    void reset(float** S) {
        X_aligned.clear();
        Y_aligned.clear();
        X_aligned.resize(0);
        Y_aligned.resize(0);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                S[i][j] = 0.0;
            }
        }

        similarity_score = 0;
        identity_score = 0;
        gap_count = 0;
    }

    // Verification of results
    bool verify(std::string file1, std::string file2) {
        std::ifstream ifs(file1);
        std::string s1(""), s2("");
        if (ifs.good()) {
            s1.assign((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>()).c_str();
        }
        ifs.close();

        ifs.open(file2);
        if (ifs.good()) {
            s2.assign((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>()).c_str();
        }
        ifs.close();

	    return (s1.size() > 0 && s1.compare(s2)==0);
    }

    // Setting up a scoring scheme without a substitution matrix
    void scoring_scheme(float match_score, float mismatch_score, float gap_penalty) {
        this->match_score = match_score;
        this->mismatch_score = mismatch_score;
        this->gap_penalty = gap_penalty;
    }
};

// Parsing arguments
void parse_args(int argc, char **argv, std::string &X, std::string &Y, std::string &output_filename, long& grain_size, int& block_size_x, int& block_size_y, int& exec_mode, bool &only_exec_times)
{
    std::string usage("Usage: --x <sequence1-filename> --y <sequence2-filename> --save-to <output-filename> --exec-mode <integer> --grain-size <integer> --block-size-x <integer> --block-size-y <integer> --print-runtime-only");

    for (int i = 0; i < argc; ++i) {
        if (std::string(argv[i]).compare("--print-runtime-only") == 0)
            only_exec_times = true;
        else if (std::string(argv[i]).compare("--x") == 0)
            X = std::string(argv[++i]);
        else if (std::string(argv[i]).compare("--exec-mode") == 0)
            exec_mode = std::stoi(argv[++i]);
        else if (std::string(argv[i]).compare("--grain-size") == 0)
            grain_size = std::stoi(argv[++i]);
        else if (std::string(argv[i]).compare("--block-size-x") == 0)
            block_size_x = std::stoi(argv[++i]);
        else if (std::string(argv[i]).compare("--block-size-y") == 0)
            block_size_y = std::stoi(argv[++i]);
        else if (std::string(argv[i]).compare("--grain-size") == 0)
            grain_size = std::stoi(argv[++i]);
        else if (std::string(argv[i]).compare("--y") == 0)
            Y = std::string(argv[++i]);
        else if (std::string(argv[i]).compare("--save-to") == 0)
            output_filename = std::string(argv[++i]);
        else if (std::string(argv[i]).compare("--help") == 0) {
            std::cout << usage << std::endl;
            exit(-1);
        }
    }
}
#endif