#include <output.h>

inline string get_id(string &name)
{
    return name.substr(0, name.find(" "));
}

void output_result_tabular(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    ofstream out(outfile);
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    for (int num_q = 0; num_q < n_query; num_q++)
    {
        vector<SWResult> &res = res_d[num_q];
        for (int t = 0; t < res.size(); t++)
        {
            out << get_id(q_names[res[t].num_q]) << "\t";
            out << get_id(res[t].s_name) << "\t";
            out << res[t].p_identity << "\t";
            out << res[t].align_length << "\t";
            out << res[t].mismatch << "\t";
            out << res[t].gap_open << "\t";
            out << res[t].begin_q << "\t";
            out << res[t].end_q << "\t";
            out << res[t].begin_s << "\t";
            out << res[t].end_s << "\t";
            out << res[t].e_value << "\t";
            out << res[t].score << "\t";
            out << "\n";
        }
    }
    out.close();
    out.clear();
}

void output_result_cast(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    ofstream out(outfile);
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    for (int num_q = 0; num_q < n_query; num_q++)
    {
        vector<SWResult> &res = res_d[num_q];
        for (int t = 0; t < res.size(); t++)
        {
            out << get_id(res[t].s_name) << "\t";
            out << get_id(q_names[res[t].num_q]) << "\t";
            out << q_names[res[t].num_q] << "\t";
            out << res[t].e_value << "\t";
            out << res[t].bitscore << "\t";
            out << res[t].score << "\t";
            out << res[t].align_length << "\t";
            out << res[t].p_identity << "\t";
            out << res[t].n_identity << "\t";
            out << res[t].mismatch << "\t";
            out << res[t].positive << "\t";
            out << res[t].gap_open << "\t";
            out << res[t].gaps << "\t";
            out << (double)res[t].positive / res[t].align_length * 100 << "\t";
            out << (double)(res[t].end_s - res[t].begin_s) / res[t].s_len * 100 << "\t";
            out << res[t].s_ori << "\t";

            out << "\n";
        }
    }
    out.close();
    out.clear();
}

void output_result_fa(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    ofstream out(outfile);
    for (int num_q = 0; num_q < n_query; num_q++)
    {
        vector<SWResult> &res = res_d[num_q];
        for (int t = 0; t < res.size(); t++)
        {
            out << ">" << res[t].s_name << "\n";
            out << res[t].s_ori << "\n";
        }
    }
    out.close();
    out.clear();
}

void output_result_a3m(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    for (int num_q = 0; num_q < n_query; num_q++)
    {
        vector<SWResult> &res = res_d[num_q];
        ofstream out(outfile + "_" + get_id(q_names[num_q]) + ".a3m");
        string q(&query[q_offsets[num_q]], &query[q_offsets[num_q + 1]] - 1);
        int qlen = q_offsets[num_q + 1] - q_offsets[num_q] - 1;
        for (int k = 0; k < qlen; k++)
        {
            q[k] += 65;
        }
        out << ">" << q_names[num_q] << "\n";
        out << q << "\n";
        for (int t = 0; t < res.size(); t++)
        {
            out << ">" << res[t].s_name << "\n";
            int qp = 0;
            while (res[t].begin_q > qp)
            {
                out << "-";
                qp++;
            }
            for (int sp = 0; sp < res[t].s.size(); sp++)
            {
                if (res[t].q[sp] == '-')
                {
                    out << (char)(res[t].s[sp] + 32);
                }
                else
                {
                    out << res[t].s[sp];
                    qp++;
                }
            }
            while (qlen > qp)
            {
                out << "-";
                qp++;
            }
            out << "\n";
        }
        out.close();
        out.clear();
    }
}

void output_result_reduce(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    ofstream out(outfile);
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    for (int num_q = 0; num_q < n_query; num_q++)
    {
        vector<SWResult> &res = res_d[num_q];
        for (int t = 0; t < res.size(); t++)
        {
            out << get_id(q_names[res[t].num_q]) << "\t";
            out << get_id(res[t].s_name) << "\t";
            out << res[t].e_value << "\t";
            out << res[t].score << "\t";
            out << res[t].bitscore << "\t";
            out << "\n";
        }
    }
    out.close();
    out.clear();
}

void output_result_align(string outfile, vector<SWResult> *res_d, const char *query, vector<uint32_t> q_offsets, vector<string> q_names)
{
    ofstream out(outfile);
    int n_query = q_names.size();
    assert(n_query == q_offsets.size() - 1);
    for (int num_q = 0; num_q < n_query; num_q++)
    {

        out << "=========================" << endl;
        string q(&query[q_offsets[num_q]], &query[q_offsets[num_q + 1]] - 1);
        int qlen = q_offsets[num_q + 1] - q_offsets[num_q] - 1;
        for (int k = 0; k < qlen; k++)
        {
            q[k] += 65;
        }
        out << endl
            << "Query " << num_q << ": " << q_names[num_q] << " #" << q << "#" << endl;

        vector<SWResult> &res = res_d[num_q];
        for (int t = 0; t < res.size(); t++)
        {
            out << endl
                << "Score = " << res[t].score << "\tE-value = " << res[t].e_value << endl;
            out << res[t].s_name << endl;

            out << "S: " << res[t].begin_s << "\t" << res[t].s << endl;
            out << "     " << "\t" << res[t].match << endl;
            out << "Q: " << res[t].begin_q << "\t" << res[t].q << endl;
        }
    }
    out.close();
    out.clear();
}