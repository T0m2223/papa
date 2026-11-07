#include <iostream>
#include <chrono>
#include <functional>
#include "helpers.hpp"
#include "implementation.hpp"

int main(int argc, char **argv)
{
    bool print_runtime_only = false;
    int exec_mode = 0; // 0. all, 1 sequential only, 2. taskloop only, 3. explicit tasks only
    long grain_size = 1; // optional parameter to use for adjusting task granularity 
    int block_size_x = 1, block_size_y = 1; // optional parameter to use for adjusting task granularity with blocks
	std::string X_filename = "X.txt", Y_filename = "Y.txt", output_filename = "aligned-sequential.txt";
	
	parse_args(argc, argv, X_filename, Y_filename, output_filename, grain_size, block_size_x, block_size_y, exec_mode, print_runtime_only);
    unsigned long entries_visited = 0, entries_visited_sequential = 0;

    SequenceInfo sinfo(X_filename, Y_filename);
    sinfo.scoring_scheme(1.0, -1.0, -2.0);
    std::cout << "Loaded X and Y sequences with sizes " << sinfo.rows -1  << " and " << sinfo.cols -1 << std::endl;
    std::cout << "Scoring scheme: match: " << sinfo.match_score  << ", mismatch " << sinfo.mismatch_score << ", gap_penalty: " << sinfo.gap_penalty << std::endl;
    std::cout << "Matrix S size: [" << sinfo.rows << "x" << sinfo.cols << "]" << std::endl;

    // allocate
    float** S = allocate(sinfo.rows,sinfo.cols, 0); // Similarity Matrix

    unsigned long expected_visited = (unsigned long)(sinfo.rows-1)*(sinfo.cols-1)+sinfo.rows+sinfo.cols-1;
    
    // sequential version
    if ( exec_mode == 1 || exec_mode < 1) {
        auto t_seq_1 = std::chrono::high_resolution_clock::now();

        entries_visited_sequential = sinfo.gpsa_sequential(S);
        
        auto t_seq_2 = std::chrono::high_resolution_clock::now();
        
        sinfo.traceback_and_save(output_filename, S);
        std::cout << "\n== Sequential version completed in " << std::chrono::duration<float>(t_seq_2 - t_seq_1).count() << " seconds." << std::endl; 
        std::cout << "   Entries visited: " << entries_visited_sequential << " " << (expected_visited ? "" : "NOT OK") << std::endl; 
        std::cout << "   Score: " << S[sinfo.rows-1][sinfo.cols-1] << ", Similarity Score: " << sinfo.similarity_score << ", Identity Score: " << sinfo.identity_score << ", Gaps: " << sinfo.gap_count << ", Length (with gaps): " << sinfo.X_aligned.size() << std::endl; 
        sinfo.reset(S);
    }
    
    // taskloop version
    if ( exec_mode == 2 || exec_mode < 1) {
        auto t_taskloop_1 = std::chrono::high_resolution_clock::now();

        entries_visited = sinfo.gpsa_taskloop(S, grain_size, block_size_x, block_size_y);
        
        auto t_taskloop_2 = std::chrono::high_resolution_clock::now();
        
        sinfo.traceback_and_save("aligned-taskloop.txt", S);

        std::cout << "\n== Taskloop version completed in " << std::chrono::duration<float>(t_taskloop_2 - t_taskloop_1).count() << " seconds." << std::endl; 
        std::cout << "   Entries visited: " << entries_visited << " " << (expected_visited == entries_visited ? "" : "NOT OK") << std::endl; 
        std::cout << "   Score: " << S[sinfo.rows-1][sinfo.cols-1] << ", Similarity Score: " << sinfo.similarity_score << ", Identity Score: " << sinfo.identity_score << ", Gaps: " << sinfo.gap_count << ", Length (with gaps): " << sinfo.X_aligned.size() << std::endl; 
        if ( exec_mode < 1 ) std::cout << "   Checking results: " << (sinfo.verify(output_filename, "aligned-taskloop.txt") ? "OK" : "NOT OK") << std::endl;
        sinfo.reset(S);
    }

    // explicit tasks versions
    if ( exec_mode == 3 || exec_mode < 1) {
        entries_visited = 0;
        auto t_tasks_1 = std::chrono::high_resolution_clock::now();

        entries_visited = sinfo.gpsa_tasks(S, grain_size, block_size_x, block_size_y);
        
        auto t_tasks_2 = std::chrono::high_resolution_clock::now();

        sinfo.traceback_and_save("aligned-tasks.txt", S);  
        
        std::cout << "\n== Explicit Tasks version completed in " << std::chrono::duration<float>(t_tasks_2 - t_tasks_1).count() << " seconds." << std::endl; 
        std::cout << "   Entries visited: " << entries_visited << " " << (entries_visited == expected_visited ? "" : "NOT OK") << std::endl; 
        std::cout << "   Score: " << S[sinfo.rows-1][sinfo.cols-1] << ", Similarity Score: " << sinfo.similarity_score << ", Identity Score: " << sinfo.identity_score << ", Gaps: " << sinfo.gap_count << ", Length (with gaps): " << sinfo.X_aligned.size() << std::endl; 
        if ( exec_mode < 1 ) std::cout << "   Checking results: " << (sinfo.verify(output_filename, "aligned-tasks.txt") ? "OK" : "NOT OK") << std::endl;
    }

    deallocate(S);

    return 0;
}