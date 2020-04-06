//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 4/4/20.
//

#include <iostream>
#include <algorithm>

#include <boost/filesystem.hpp>

#include <gflags/gflags.h>

#include <sdsl/construct.hpp>

DEFINE_string(text, "", "Text file. (MANDATORY)");
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

using namespace sdsl;
using namespace std;


const char* KEY_BWT_RUNS_FIRST = "bwt_runs_first";
const char* KEY_BWT_RUNS_LAST = "bwt_runs_last";
const char* KEY_TEXT_BWT_RUNS_FIRST = "text_bwt_runs_first";
const char* KEY_TEXT_BWT_RUNS_LAST = "text_bwt_runs_last";


//! Constructs the Burrows and Wheeler Transform (BWT) from text over byte- or integer-alphabet and suffix array.
/*!	The algorithm constructs the BWT and stores it to disk.
 *  \tparam t_width Width of the text. 0==integer alphabet, 8=byte alphabet.
 *  \param config	Reference to cache configuration
 *  \par Space complexity
 *		\f$ n \log \sigma \f$ bits
 *  \pre Text and Suffix array exist in the cache. Keys:
 *         * conf::KEY_TEXT for t_width=8 or conf::KEY_TEXT_INT for t_width=0
 *         * conf::KEY_SA
 *  \post BWT exist in the cache. Key
 *         * conf::KEY_BWT for t_width=8 or conf::KEY_BWT_INT for t_width=0
 */
template<uint8_t t_width>
void construct_bwt_and_runs(cache_config& config)
{
    static_assert(t_width == 0 or t_width == 8 , "construct_bwt_and_runs: width must be `0` for integer alphabet and `8` for byte alphabet");

    typedef int_vector<>::size_type size_type;
    typedef int_vector<t_width> text_type;
    typedef int_vector_buffer<t_width> bwt_type;
    const char* KEY_TEXT = key_text_trait<t_width>::KEY_TEXT;
    const char* KEY_BWT = key_bwt_trait<t_width>::KEY_BWT;

    //  (1) Load text from disk
    text_type text;
    load_from_cache(text, KEY_TEXT, config);
    size_type n = text.size();
    uint8_t bwt_width = text.width();

    //  (2) Prepare to stream SA from disc and BWT to disc
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
    int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, config), std::ios::in, buffer_size);
    std::string bwt_file = cache_file_name(KEY_BWT, config);
    bwt_type bwt_buf(bwt_file, std::ios::out, buffer_size, bwt_width);

    int_vector_buffer<> bwt_runs_first_buf(cache_file_name(KEY_BWT_RUNS_FIRST, config), std::ios::out);
    int_vector_buffer<> text_bwt_runs_first_buf(cache_file_name(KEY_TEXT_BWT_RUNS_FIRST, config), std::ios::out);
    int_vector_buffer<> text_bwt_runs_last_buf(cache_file_name(KEY_TEXT_BWT_RUNS_LAST, config), std::ios::out);

    //  (3) Construct BWT sequentially by streaming SA and random access to text
    size_type to_add[2] = {(size_type) -1, n - 1};
    auto get_text_position = [&sa_buf, &to_add](auto i) {
        auto pos = sa_buf[i];
        return pos + to_add[pos == 0];
    };

    // First BWT value
    bwt_buf[0] = text[ get_text_position(0) ];
    uint8_t prev_c_bwt = bwt_buf[0];

    size_t nruns = 0; // # of BWT runs

    // First position starts the first BWT run.
    bwt_runs_first_buf[nruns] = 0;
    text_bwt_runs_first_buf[nruns] = get_text_position(0);

    auto last_pos = n - 1;
    for (size_type i = 1; i < last_pos; ++i) {
        bwt_buf[i] = text[ get_text_position(i) ];

        uint8_t c_bwt = bwt_buf[i];
        if (prev_c_bwt != c_bwt) {
            // Last position of the current BWT run
            text_bwt_runs_last_buf[nruns] = get_text_position(i - 1);

            ++nruns;

            // First position of the next BWT run
            bwt_runs_first_buf[nruns] = i;
            text_bwt_runs_first_buf[nruns] = get_text_position(i);

            prev_c_bwt = c_bwt;
        }
    }

    // Last BWT value
    bwt_buf[last_pos] = text[ get_text_position(last_pos) ];

    // Last position ends the last BWT run
    text_bwt_runs_last_buf[nruns] = get_text_position(last_pos);

    cout << "# BWT runs = " << nruns << endl;

    text_bwt_runs_first_buf.close();
    text_bwt_runs_last_buf.close();
    bwt_runs_first_buf.close();

    bwt_buf.close();

    register_cache_file(KEY_BWT, config);
}

int main(int argc, char **argv) {
    gflags::SetUsageMessage("This program calculates the SA and BWT for the given text.");
    gflags::AllowCommandLineReparsing();
    gflags::ParseCommandLineFlags(&argc, &argv, false);

    if (FLAGS_text.empty()) {
        std::cerr << "Command-line error!!!" << std::endl;
        return 1;
    }

    string data_path = FLAGS_text;

    cache_config config(false, ".", util::basename(FLAGS_text));

    if (!cache_file_exists(conf::KEY_TEXT, config)) {
        int_vector<8> text;
        {
            std::string input;
            {
                std::ifstream fs(data_path);
                std::stringstream buffer;
                buffer << fs.rdbuf();

                input = buffer.str();
            }

            cout << "Text length = " << input.size() << endl;
//        if (any_of(text.begin(), text.end(), [](auto symbol) { return symbol == 0 || symbol == 1; })) {
//            std::cout << "ERROR: Input data contains reserved characters {0x0, 0x1}";
//            exit(1);
//        }

            text.resize(input.size() + 1);

            replace_copy(input.begin(), input.end(), text.begin(), 0, 2);

//        append_zero_symbol(text);
            text[text.size() - 1] = 0;
        }

        store_to_cache(text, conf::KEY_TEXT, config);
        util::clear(text);
    }

    if (!cache_file_exists(conf::KEY_SA, config)) {
        cout << "Calculate Suffix Array ... " << endl;
        construct_config::byte_algo_sa = FLAGS_sais ? SE_SAIS
                                                    : LIBDIVSUFSORT; // or LIBDIVSUFSORT for less space-efficient but faster construction
        construct_sa<8>(config);
        cout << "DONE" << endl;
    }

    if (!cache_file_exists(conf::KEY_BWT, config)) {
        cout << "Calculate BWT ... " << endl;
//        construct_bwt<8>(config);
        construct_bwt_and_runs<8>(config);
        cout << "DONE" << endl;
    }

    return 0;
}
