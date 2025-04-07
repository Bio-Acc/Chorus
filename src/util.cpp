#include "util.h"

#define MASK5(v) (v & 0b11111)

ArgumentParser arg_parser;
int seed_length;
int qit_length;
int qit_width;
int max_display_align;
int band_width;
int num_threads;
int min_score;
double max_evalue;
bool detailed_alignment = false;
int gpu_device;
vector<Condition> must_include;
bool has_must_include = false;

// StripedSmithWaterman::Aligner aligner((-1) * SCORE_GAP, (-1) * SCORE_GAP_EXT);
// StripedSmithWaterman::Filter filter;

double timeuse(struct timeval start_time, struct timeval end_time)
{
    return (end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_usec - start_time.tv_usec) / 1000000.0;
}

int get_pos(size_t idx, const size_t *offsets, size_t n_subj)
{
    auto p = upper_bound(offsets, offsets + n_subj, idx) - 1;
    return p - offsets;
}

uint32_t load_fasta(const char *file, char *&str, vector<uint32_t> &offsets, vector<string> &names)
{

    std::ifstream input(file);
    if (!input.good())
    {
        std::cerr << "Error opening '" << file << "'. Bailing out." << std::endl;
        return 0;
    }
    str = nullptr;
    size_t ptr = 0;

    std::string line, name, content;
    while (std::getline(input, line).good())
    {
        if (line.empty() || line[0] == '>')
        { // Identifier marker
            if (!name.empty())
            { // Print out what we read from the last entry
                // std::cout << name << " : " << content.size() << std::endl;
                str = (char *)realloc(str, ptr + content.size() + 1);
                for (int k = 0; k < content.size(); k++)
                {
                    if (content[k] >= 65 && content[k] <= 90)
                        content[k] -= 65;
                    else if (content[k] >= 97 && content[k] <= 122)
                        content[k] -= 97;
                    else
                    {
                        cout << "warning:" << content[k] << " " << (int)content[k] << endl;
                        content[k] = ILLEGAL_WORD;
                    }
                }
                content.copy(str + ptr, content.size(), 0);
                offsets.emplace_back(ptr);
                names.emplace_back(name);
                ptr += content.size();
                str[ptr++] = END_SIGNAL;
                name.clear();
            }
            if (!line.empty())
            {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty())
        {
            if (line.find(' ') != std::string::npos)
            { // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            }
            else
            {
                content += line;
            }
        }
    }
    if (!name.empty())
    { // Print out what we read from the last entry
        // std::cout << name << " : " << content << std::endl;
        content +=line;
        str = (char *)realloc(str, ptr + content.size() + 1);
        for (int k = 0; k < content.size(); k++)
        {
            if (content[k] >= 65 && content[k] <= 90)
                content[k] -= 65;
            else
            {
                cout << "warning:" << content[k] << " " << (int)content[k] << endl;
                content[k] = ILLEGAL_WORD;
            }
        }
        content.copy(str + ptr, content.size(), 0);
        offsets.emplace_back(ptr);
        names.emplace_back(name);
        ptr += content.size();
        str[ptr++] = END_SIGNAL;
    }

    cout << "Load " << file << " successful. Num of seq = " << names.size() << ". Total length = " << ptr << endl;
    input.close();
    input.clear();
    return ptr;
}

size_t load_offsets(string f_name, size_t *&s_offsets, size_t *&n_offsets)
{
    int fds = open((f_name + ".sofs").data(), O_RDONLY);
    if (fds == -1)
    {
        std::cerr << "Error opening '" << f_name << "'.sofs Bailing out." << std::endl;
        return 0;
    }
    size_t len = lseek(fds, 0, SEEK_END);
    size_t *ss_offsets = (size_t *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fds, 0);
    close(fds);

    s_offsets = (size_t *)malloc(len);
    memcpy(s_offsets, ss_offsets, len);
    munmap(ss_offsets, len);

    int fdn = open((f_name + ".nofs").data(), O_RDONLY);
    if (fdn == -1)
    {
        std::cerr << "Error opening '" << f_name << "'.nofs Bailing out." << std::endl;
        return 0;
    }
    size_t *nn_offsets = (size_t *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fdn, 0);
    close(fdn);

    n_offsets = (size_t *)malloc(len);
    memcpy(n_offsets, nn_offsets, len);
    munmap(nn_offsets, len);
    return len;
}

// size_t load_name(string f_name, int num, char *&s_names)
// {
//     int fd = open((f_name + "_" + to_string(num) + ".name").data(), O_RDONLY);
//     if (fd == -1)
//     {
//         std::cerr << "Error opening '" << f_name << "'.name Bailing out." << std::endl;
//         return 0;
//     }
//     size_t len = lseek(fd, 0, SEEK_END);
//     char *ss_names = (char *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
//     close(fd);

//     s_names = (char *)malloc(len);
//     memcpy(s_names, ss_names, len);
//     munmap(ss_names, len);
//     return len;
// }

int check_db(const char *db, size_t &max_size, size_t &total_size)
{
    int i = 0;
    while (true)
    {
        string f_name(db);
        size_t this_size = 0;
        for (int s = 0; s < NUM_STREAM; s++)
        {
            std::ifstream input(f_name + "_" + to_string(i) + "_" + to_string(s) + ".seq");
            if (!input.good())
            {
                cout << "Found db " << f_name << ", total parts = " << i << endl;
                return i;
            }
            input.seekg(0, ios::end);
            size_t len = input.tellg();
            this_size += len - 20;
            input.close();
            input.clear();
        }
        total_size += this_size;
        max_size = max_size < this_size ? this_size : max_size;
        i++;
    }
    return 0;
}

// string get_name(string db_name, size_t offset, size_t len)
// {
//     // db_name.clear();
//     // db_name.seekg(offset, ios::beg);
//     // string name;
//     // std::getline(db_name, name);
//     // assert(name.length() == len);
//     int fd = open(db_name.data(), O_RDONLY);
//     if (fd == -1)
//     {
//         std::cerr << "Error opening '" << db_name << "'. Bailing out." << std::endl;
//         return "";
//     }
//     size_t file_len = offset + len;
//     char *str1 = (char *)mmap(NULL, file_len, PROT_READ, MAP_PRIVATE, fd, 0);
//     if (str1 == MAP_FAILED)
//         cout << "map failed" << endl;
//     close(fd);
//     string name(str1 + offset, len);
//     munmap(str1, file_len);
//     // cout<<name<<endl;;
//     return name;
// }

void load_seq(string db, int num, char *&str, size_t &len)
{
    string f_name(db);
    f_name = f_name + "_" + to_string(num) + ".seq";

    int fd = open(f_name.data(), O_RDONLY);
    if (fd == -1)
    {
        std::cerr << "Error opening '" << f_name << "'. Bailing out." << std::endl;
        return;
    }

    len = lseek(fd, 0, SEEK_END);
    // 建立内存映射
    char *str1 = (char *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    // str=(char*)malloc(len);
    // cudaHostAlloc((void **)&str, len, cudaHostAllocDefault);
    // memcpy(str, str1, len);
    // munmap(str1, len);
    str = str1;
    // cout << "Load db " << db << num << " successful. Total length = " << len << endl;
}

void proceed_result(vector<SWResult> *res_d, vector<SWResult> &res_t, const char *query, const char *subj, QueryGroup &q_group, const char *s_name, const size_t *s_offsets, const size_t *sn_offsets, const size_t s_num, size_t total_db_size)
{
    set<int> modified_resq;

    // struct timeval t_start, t_end;
    // gettimeofday(&t_start, NULL);
    // size_t s_num = load_offsets(db_name + "_" + to_string(db_num) + "_" + to_string(stream_num), s_offsets, sn_offsets);

    // sort(res_t.begin(), res_t.end(), [&](const SWResult &sw1, const SWResult &sw2)
    // {
    //     if (sw1.num_q==sw2.num_q)
    //         return (sw1.score>sw2.score);
    //     else return sw1.num_q<sw2.num_q; });

    size_t n_subj = s_num / sizeof(size_t);

    // vector<SWResult> res_list;
    // res_list.resize(res_t.size());

#pragma omp parallel for
    for (int t = 0; t < res_t.size(); t++)
    {
        if (!res_t[t].report)
            continue;
        // double score = (E_lambda * res_t[t].score - log(E_k)) / (0.69314718055995);
        double score = res_t[t].score;
        if (score >= min_score)
        {
            // if (res_t[t].s_res.size()<=0)
            // {
            //     cout << score << endl;
            // }
            // cout << "s_begin" << t<<":"<<endl;
            // cout <<res_t[t].s_res[0]<<endl;
            int num_s = get_pos(res_t[t].s_res[0], s_offsets, n_subj);
            int num_q = res_t[t].num_q;
            res_t[t].num_q = q_group.id[num_q];
            uint32_t s_len = s_offsets[num_s + 1] - s_offsets[num_s] - 1;

            double e_value = total_db_size * (q_group.length[num_q]) * pow(2, -res_t[t].bitscore);
            if (e_value <= max_evalue)
            {
                res_t[t].e_value = e_value;
                res_t[t].sn_offset = sn_offsets[num_s];
                res_t[t].sn_len = sn_offsets[num_s + 1] - sn_offsets[num_s] - 1;

                string sn(s_name + res_t[t].sn_offset, res_t[t].sn_len);
                res_t[t].s_name = sn;

                res_t[t].begin_s = res_t[t].s_res[0] - s_offsets[num_s];
                res_t[t].begin_q = res_t[t].q_res[0] - q_group.offset[num_q];

                size_t len = res_t[t].align_length;
                res_t[t].end_q = res_t[t].q_res[len - 1] - q_group.offset[num_q];
                res_t[t].end_s = res_t[t].s_res[len - 1] - s_offsets[num_s];

                res_t[t].s_len = s_len;
            }
            else
            {
                res_t[t].report=false;
                res_t[t].q.clear();
                res_t[t].s.clear();
                res_t[t].s_ori.clear();
                res_t[t].match.clear();
            }
        }
        else
        {
            res_t[t].report=false;
            res_t[t].q.clear();
            res_t[t].s.clear();
            res_t[t].s_ori.clear();
            res_t[t].match.clear();
        }
        res_t[t].s_res.clear();
        res_t[t].q_res.clear();

        // res_t[t].s_res.shrink_to_fit();
        // res_t[t].q_res.shrink_to_fit();
        
    }
    for (int t = 0; t < res_t.size(); t++)
    {
        if (!res_t[t].report)
            continue;
        // SWResult& res = res_t[t];
        res_d[res_t[t].num_q].emplace_back(res_t[t]);
        push_heap(res_d[res_t[t].num_q].begin(), res_d[res_t[t].num_q].end(), [&](const SWResult &sw1, const SWResult &sw2)
                          { return (sw1.e_value == sw2.e_value) ? (sw1.score > sw2.score) : (sw1.e_value < sw2.e_value); });
                // res_tmp.emplace_back(res);
        modified_resq.insert(res_t[t].num_q);
    }
    // munmap(s_offsets, s_num);
    // munmap(sn_offsets, s_num);
    // free(s_offsets);
    // free(sn_offsets);
    // free(s_names);

    // int max_display_align;
    // get_arg("max-output-align", max_display_align, 0);
    // cout << "res size ";

    // sort(res_tmp.begin(), res_tmp.end(), [&](const DisplayResult &sw1, const DisplayResult &sw2)
    //          { return (sw1.num_q == sw2.num_q) ? (sw1.e_value < sw2.e_value):(sw1.num_q < sw2.num_q); });

    // if (!modified_resq.empty()) {
    //     for (auto i = modified_resq.begin(); i != modified_resq.end(); ++i) {
    //         sort(i->second.begin(), i->second.end(), [&](const DisplayResult &sw1, const DisplayResult &sw2)
    //          { return (sw1.e_value == sw2.e_value) ? (sw1.score > sw2.score) : (sw1.e_value < sw2.e_value); });
    //         int k = res_d[i->first].size();
    //         res_d[i->first].resize(k + i->second.size());
    //         merge(res_d[i->first], k, i->second, i->second.size());
    //         if (max_display_align > 0 && res_d[i->first].size() > max_display_align)
    //         {
    //             res_d[i->first].resize(max_display_align);
    //             assert(res_d[i->first].size() == max_display_align);
    //         }
    //     }
    // }
    vector<int> modified_resq_list;
    modified_resq_list.assign(modified_resq.begin(), modified_resq.end());//用assign实现set转vector

#pragma omp parallel for
    for (int i=0;i<modified_resq_list.size();i++)
    {
        int mq = modified_resq_list[i];
        vector<SWResult> &res = res_d[mq];

        // sort(res.begin(), res.end(), [&](const DisplayResult &sw1, const DisplayResult &sw2)
        //      { return (sw1.e_value == sw2.e_value) ? (sw1.score > sw2.score) : (sw1.e_value < sw2.e_value); });
        if (max_display_align > 0 && res.size() > max_display_align)
        {
            int pop_size = res.size() - max_display_align;
            for (int i = 0; i < pop_size; i++)
            {
                pop_heap(res.begin(), res.end() - i, [&](const SWResult &sw1, const SWResult &sw2)
                         { return (sw1.e_value == sw2.e_value) ? (sw1.score > sw2.score) : (sw1.e_value < sw2.e_value); });
                // res.pop_back();
            }
            // res.erase(res.begin() + max_display_align, res.end());
            res.resize(max_display_align);
            assert(res.size() == max_display_align);
        }
        // cout <<"("<<mq<<"," <<res.size() <<")";
        // if (max_display_align > 0 && res.size() > max_display_align)
        // {
        //     std::cout << "error" << res.size() <<" "<<mq<<std::endl;
        //     exit(-1);
        // }
    }
    // cout<<endl;
}

inline char get_char(const char *s, size_t offset)
{
    offset *= 5;
    return MASK5((unsigned)((*((uint16_t *)&(s[offset >> 3]))) >> (offset & 7)));
}

void get_arg(const char *name, int &v, int default_v)
{
    string arg = arg_parser.retrieve<string>(name);
    if (arg.compare("") == 0)
        v = default_v;
    else
        v = std::stoi(arg);
}

void get_arg(const char *name, double &v, double default_v)
{
    string arg = arg_parser.retrieve<string>(name);
    if (arg.compare("") == 0)
        v = default_v;
    else
        v = std::stod(arg);
}

void get_arg(const char *name, string &v)
{
    v = arg_parser.retrieve<string>(name);
}

void get_arg(const char *name, vector<string> &v)
{
    v = arg_parser.retrieve<vector<string>>(name);
}

void load_must_include()
{
    has_must_include = arg_parser.exists("must-include");
    if (has_must_include)
    {
        vector<string> argvs;
        get_arg("must-include", argvs);
        int n = argvs.size();
        if (n==0)
        {
            has_must_include = false;
            return;
        }
        if (n % 2 != 0)
        {
            cerr << "Error param: --must-incluce [regular expression] [match times]" << endl;
            exit(-1);
        }
        n = n / 2;
        for (int i = 0; i < n; i++)
        {
            Condition c;
            c.re = argvs[i * 2];
            c.t = stoi(argvs[i * 2 + 1]);
            must_include.push_back(c);
        }
        // string fname;
        // get_arg("must-include", fname)
        // std::ifstream input(fname);
        // if (!input.good())
        // {
        //     std::cerr << "Error opening '" << file << "'. Bailing out." << std::endl;
        //     exit(-1);
        // }
        // std::string line;
        // while (std::getline(input, line).good())
        // {

        // }
    }
}

bool check_include(const char *s, size_t begin_p, size_t end_p)
{
    char str[end_p - begin_p + 1] = {0};
    for (size_t i = begin_p; i < end_p; i++)
    {
        str[i - begin_p] = get_char(s, i) + 65;
    }
    string ss(str, end_p - begin_p);
    return check_include(ss);
}

bool check_include(string str)
{
    for (auto c : must_include)
    {
        string subj(str);
        smatch m;
        regex e(c.re);
        int times = 0;
        while (regex_search(subj, m, e))
        {
            times++;
            if (times >= c.t)
                break;
            subj = m.suffix().str();
        }
        if (times < c.t)
            return false;
    }
    return true;
}