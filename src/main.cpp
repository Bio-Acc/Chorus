
#include "blastp.h"

void show_help_info(int argc, const char **argv)
{
    if (argc <= 1 || !strcmp(argv[1], "--help"))
    {
        cout << "usage: " << argv[0] << " -q <query.fasta> -d <db1> [db2 ... (optional)] -o <output file> [options]" << endl;
        cout << endl;
        cout << "Algorithm options:" << endl;
        cout << endl;
        cout << "\t"
             << "--outfmt"
             << "\t\t"
             << "INT"
             << "\t"
             << "Output format: 0: m8 tabular (like blast+ outfmt 6); 1: detailed alignment; 2: tabular and full ref seqs (fasta); 3: Opfi format; 4: a3m format; 5: Diamond benchmark format. [" << D_OUTFMT << "]" << endl;
        cout << "\t"
             << "-l, --filter-level"
             << "\t"
             << "INT"
             << "\t"
             << "The filter level to pass a HSP candidate to local alignment. The smaller the value, the more sensitive the result. [" << D_FILTER_LEVEL << "]" << endl;
        cout << "\t"
             << "--min-score"
             << "\t\t"
             << "INT"
             << "\t"
             << "The minimum score for each alignment to display. [" << D_MIN_SCORE << "]" << endl;
        cout << "\t"
             << "-e, --max-evalue"
             << "\t"
             << "FLOAT"
             << "\t"
             << "The maximum expectation value for each alignment to display. [" << D_MAX_EVALUE << "]" << endl;
        cout << "\t"
             << "--max-output-align"
             << "\t"
             << "INT"
             << "\t"
             << "The maximum number of alignments to display for each query (0 means no limit). [" << D_MAX_DISPLAY_ALIGN << "]" << endl;
        cout << "\t"
             << "--num-threads"
             << "\t\t"
             << "INT"
             << "\t"
             << "Number of threads for CPU local alignment. [" << thread::hardware_concurrency() << "]" << endl;
        cout << "\t"
             << "-k, --seed-length"
             << "\t"
             << "INT"
             << "\t"
             << "Length of k-mer (3,4,5). [" << D_SEED_LENGTH << "]" << endl;
        cout << "\t"
             << "-w, --qit-width"
             << "\t\t"
             << "INT"
             << "\t"
             << "Width of query k-mer index table. [" << D_QIT_WIDTH << "]" << endl;
        cout << "\t"
             << "-h, --hash-size"
             << "\t\t"
             << "INT"
             << "\t"
             << "Voting hash table size ratio. The larger the value, the larger the GPU voting hash table (doubly for each increment). (0: maybe OK; 1: OK for most query seqs; 2: redundancy) [" << D_HASH_RATIO << "]" << endl;
        cout << "\t"
             << "--band-width"
             << "\t\t"
             << "INT"
             << "\t"
             << "The bandwidth in banded local alignment. [" << D_BAND_WIDTH_CPU << "]" << endl;
        cout << "\t"
             << "--must-include"
             << "\t\t"
             << "STRING INT"
             << "\t"
             << "Requires the ref sequence must match a regular expression for N times. (e.g. --must-include \"R[A-Z]{4,7}H\" 2)" << endl;
        cout << endl;
        exit(0);
    }
}

void show_args()
{
    string argv_query, argv_out;
    vector<string> argv_dbs;
    get_arg("query", argv_query);
    get_arg("db", argv_dbs);
    get_arg("out", argv_out);

    cout << "Alignment Type" << "\t";
#ifdef USE_GPU_SW
    cout << "GPU SW" << endl;
#elif defined(GLF_GPU_SW)
     cout << "GLF SW" << endl;
#else
    cout << "CPU BANDED SW" << endl;
#endif

    cout << "Query"
         << "\t" << argv_query << endl;
    cout << "Ref DB" << endl;
    for (int i = 0; i < argv_dbs.size(); i++)
    {
        cout << "\t" << argv_dbs[i] << endl;
    }
    cout << "Output"
         << "\t" << argv_out << endl;

    int filter_level, hashtable_size_ratio;

    get_arg("filter-level", filter_level, D_FILTER_LEVEL);
    get_arg("hash-size", hashtable_size_ratio, D_HASH_RATIO);
    load_must_include();

    cout << "K-mer Length"
         << "\t" << seed_length << endl;
    cout << "Query Index Table Width"
         << "\t" << qit_width << endl;
    cout << "Max Output Align"
         << "\t" << max_display_align << endl;
    cout << "Min Score"
         << "\t" << min_score << endl;
    cout << "Max E-value"
         << "\t" << max_evalue << endl;
    cout << "Filter Level"
         << "\t" << filter_level << endl;

    cout << "Num Threads"
         << "\t" << num_threads << endl;
    cout << "Num Streams"
         << "\t" << NUM_STREAM << endl;
    cout << "Voting Hashtable Size Ratio"
         << "\t" << hashtable_size_ratio << endl;
    cout << "SW Band Width"
         << "\t" << band_width << endl;

    cout << "Scoring Matrix"
         << "\t"
         << "BLOSUM62" << endl;
    cout << "K in E_value"
         << "\t" << E_k << endl;
    cout << "Lambda in E_value"
         << "\t" << E_lambda << endl;

    if (has_must_include)
    {
        cout << "Must include:" << endl;
        for (auto c : must_include)
        {
            cout << "\t" << c.re << "\t" << c.t << endl;
        }
    }
}

int main(int argc, const char **argv)
{
    // ios::sync_with_stdio(false);
    show_help_info(argc, argv);
    arg_parser.addArgument("-d", "--db", '+', false);
    arg_parser.addArgument("-q", "--query", 1, false);
    arg_parser.addArgument("-o", "--out", 1, false);
    arg_parser.addArgument("--max-output-align", 1, true);
    arg_parser.addArgument("--min-score", 1, true);
    arg_parser.addArgument("-e", "--max-evalue", 1, true);
    arg_parser.addArgument("--num-threads", 1, true);
    arg_parser.addArgument("--band-width", 1, true);
    arg_parser.addArgument("-l", "--filter-level", 1, true);
    arg_parser.addArgument("-k", "--seed-length", 1, true);
    arg_parser.addArgument("-w", "--qit-width", 1, true);
    arg_parser.addArgument("-h", "--hash-size", 1, true); // 0: maybe ok, 1: ok, 2: too large
    arg_parser.addArgument("--outfmt", 1, true);
    arg_parser.addArgument("--must-include", '+', true);
    arg_parser.parse(argc, argv);

    string argv_query, argv_out;
    vector<string> argv_dbs;
    get_arg("query", argv_query);
    get_arg("db", argv_dbs);
    get_arg("out", argv_out);
    get_arg("num-threads", num_threads, thread::hardware_concurrency());
    get_arg("seed-length", seed_length, D_SEED_LENGTH);
    get_arg("qit-width", qit_width, D_QIT_WIDTH);
    qit_length = (int)pow(32, seed_length);
    get_arg("max-output-align", max_display_align, D_MAX_DISPLAY_ALIGN);
    get_arg("min-score", min_score, D_MIN_SCORE);
    get_arg("max-evalue", max_evalue, D_MAX_EVALUE);
#ifdef USE_GPU_SW
     get_arg("band-width", band_width, D_BAND_WIDTH_GPU);
#else
     get_arg("band-width", band_width, D_BAND_WIDTH_CPU);
#endif
    show_args();

    int outfmt;
    get_arg("outfmt", outfmt, D_OUTFMT);
    switch (outfmt)
    {
    case 0:
        break;
    case 1:
    case 2:
    case 3:
    case 4:
        detailed_alignment = true;
        break;
    case 5:
        break;
    default:
        break;
    }

    omp_set_num_threads(num_threads);

    struct timeval tt_start, tt_end;
    gettimeofday(&tt_start, NULL);
    blastp(argv_query, argv_dbs, argv_out);
    gettimeofday(&tt_end, NULL);
    cout << "Total Time:\t" << timeuse(tt_start, tt_end) << endl;

    return 0;
}