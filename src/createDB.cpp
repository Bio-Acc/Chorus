
#include "params.h"
#include "util.h"


#define MAX_QUERY_LENGTH 10000

#define BUFFER_SIZE 20

using namespace std;

uint32_t count_3mer[32][32][32]={0};
uint32_t count_4mer[32][32][32][32]={0};
uint32_t count_5mer[32][32][32][32][32]={0};

void help(int argc, const char **argv)
{
    if (argc <= 1 || !strcmp(argv[1], "--help"))
    {
        cout << "usage: " << argv[0] << " <protein_db.fasta> <output db> <batch size (GB)>" << endl;
        exit(0);
    }
}

// void visualize_alphabet(bool alphabet[128])
// {
//     for (int i = 0; i < 128; i++)
//     {
//         if (alphabet[i])
//             printf("%d:%c\n", i, i);
//     }
// }

// inline char get_char(char *s, size_t offset)
// {
//     size_t n_bit = offset * 5;
//     return MASK5((unsigned)((*((uint16_t *)&(s[n_bit >> 3]))) >> (n_bit & 7))) + 65;
// }

inline void buffer_insert(char *buffer, int &buffer_bit, char ch, int& content_gaps)
{
    if (ch >= 65 && ch <= 90)
        ch -= 65;
    else if (ch >= 97 && ch <= 122)
        ch -= 97;
    else if (ch == '-')
    {
        content_gaps++;
        return;
    }
    else if (ch != -1)
        ch = -2;
    unsigned buffer_byte = buffer_bit >> 3;
    unsigned buffer_byte_bit = buffer_bit & 7;
    buffer[buffer_byte] |= (MASK5(ch) << (buffer_byte_bit));
    if (buffer_byte_bit > 3)
    {
        buffer[buffer_byte + 1] |= (MASK5(ch) >> (8 - buffer_byte_bit));
    }
    buffer_bit += 5;
}

// 0-25: A-Z, -1: end of a sequence, -23: illegal character

void create_db(const char *fa_file, const char *out_path, size_t max_size)
{
    std::ifstream input(fa_file);
    if (!input.good())
    {
        std::cerr << "Error opening '" << fa_file << "'. Bailing out." << std::endl;
        return;
    }

    input.seekg(0, ios::end);
    size_t f_size = input.tellg();
    input.seekg(0);

    std::string line, name, content;
    size_t ptr = 0;
    int content_gaps = 0;
    int buffer_bit = 0;
    char buffer[BUFFER_SIZE] = {0};
    size_t n_ptr = 0;
    size_t db_i = 0;
    size_t stream_i = 0;
    size_t s_num = 0;
    size_t db_pointer = input.tellg();
    size_t db_begin = input.tellg();

    string f_name = out_path;
    // string db_name = f_name + "_infos.sqlite";
    ofstream *db_seq = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".seq");
    ofstream *db_name = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".name");
    ofstream *db_sofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".sofs");
    ofstream *db_nofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".nofs");

    // sqlite3 *db_info = create_db(db_name.data());

    // sqlite3_exec(db_info, "begin;", 0, 0, 0);
    while (std::getline(input, line).good())
    {
        if (line.empty() || line[0] == '>' || line[0] == '#')
        {
            if (!name.empty())
            {
                if (ptr + content.size() >= max_size)
                {
                    // *db_seq << buffer;
                    while (buffer_bit < BUFFER_SIZE * 8)
                    {
                        buffer_insert(buffer, buffer_bit, -1, content_gaps);
                    }
                    db_seq->write(buffer, BUFFER_SIZE);
                    memset(buffer, 0xFF, BUFFER_SIZE);
                    db_seq->write(buffer, BUFFER_SIZE);
                    buffer_bit = 0;
                    db_sofs->write((char *)&ptr, sizeof(size_t));
                    db_nofs->write((char *)&n_ptr, sizeof(size_t));
                    db_seq->close();
                    db_seq->clear();
                    delete db_seq;
                    db_name->close();
                    db_name->clear();
                    delete db_name;
                    // sqlite3_exec(db_info, "commit;", 0, 0, 0);
                    db_sofs->close();
                    db_sofs->clear();
                    delete db_sofs;
                    db_nofs->close();
                    db_nofs->clear();
                    delete db_nofs;
                    cout << (int)((double)input.tellg() / f_size * 100) << "%\t Create db " << db_i << "_" << stream_i << " successful! Num. sequence = " << s_num << ", Sequence length = " << ptr << endl;
                    stream_i++;
                    if (stream_i >= NUM_STREAM)
                    {
                        db_i++;
                        stream_i = 0;
                        db_begin = db_pointer - 1;
                    }
                    db_seq = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".seq");
                    db_name = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".name");
                    db_sofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".sofs");
                    db_nofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".nofs");
                    // sqlite3_exec(db_info, "begin;", 0, 0, 0);
                    ptr = 0;
                    n_ptr = 0;
                    s_num = 0;
                    memset(buffer, 0, BUFFER_SIZE);
                }
                content.push_back(-1);
                content_gaps = 0;
                for (int j = 0; j < content.size(); j++)
                {
                    char ch = content[j];
                    buffer_insert(buffer, buffer_bit, ch, content_gaps);
                    if (buffer_bit >= BUFFER_SIZE * 8)
                    {
                        assert(buffer_bit == BUFFER_SIZE * 8);
                        // *db_seq << buffer;
                        db_seq->write(buffer, BUFFER_SIZE);
                        db_seq->clear();
                        memset(buffer, 0, BUFFER_SIZE);
                        buffer_bit = 0;
                    }
                }
                *db_name << name << '\n';
                // insert_seq(db_info, db_i, stream_i, n_ptr, name);
                db_sofs->write((char *)&ptr, sizeof(size_t));
                db_nofs->write((char *)&n_ptr, sizeof(size_t));
                s_num++;
                // *db_sofs << ptr << endl;
                // *db_nofs << n_ptr << endl;
                ptr += content.size()- content_gaps;
                n_ptr += name.size() + 1;
                name.clear();
            }
            if (!line.empty())
            {
                if (line[0] == '>')
                {
                    db_pointer = input.tellg();
                    db_pointer -= line.length();
                    name = line.substr(1);
                }
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
        if (ptr + content.size() >= max_size)
        {
            // *db_seq << buffer;
            while (buffer_bit < BUFFER_SIZE * 8)
            {
                buffer_insert(buffer, buffer_bit, -1, content_gaps);
            }
            db_seq->write(buffer, BUFFER_SIZE);
            memset(buffer, 0xFF, BUFFER_SIZE);
            db_seq->write(buffer, BUFFER_SIZE);
            buffer_bit = 0;
            db_sofs->write((char *)&ptr, sizeof(size_t));
            db_nofs->write((char *)&n_ptr, sizeof(size_t));
            db_seq->close();
            db_seq->clear();
            delete db_seq;
            db_name->close();
            db_name->clear();
            delete db_name;
            // sqlite3_exec(db_info, "commit;", 0, 0, 0);
            db_sofs->close();
            db_sofs->clear();
            delete db_sofs;
            db_nofs->close();
            db_nofs->clear();
            delete db_nofs;
            cout << (int)((double)input.tellg() / f_size * 100) << "%\t Create db " << db_i << "_" << stream_i << " successful! Num. sequence = " << s_num << ", Sequence length = " << ptr << endl;
            stream_i++;
            if (stream_i >= NUM_STREAM)
            {
                db_i++;
                stream_i = 0;
                db_begin = db_pointer - 1;
            }
            db_seq = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".seq");
            db_name = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".name");
            db_sofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".sofs");
            db_nofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".nofs");
            // sqlite3_exec(db_info, "begin;", 0, 0, 0);
            // cout << name << endl;
            ptr = 0;
            n_ptr = 0;
            s_num = 0;
            memset(buffer, 0, BUFFER_SIZE);
        }
        content.push_back(-1);
        for (int j = 0; j < content.size(); j++)
        {
            char ch = content[j];

            if (buffer_bit >= BUFFER_SIZE * 8)
            {
                assert(buffer_bit == BUFFER_SIZE * 8);
                // *db_seq << buffer;
                db_seq->write(buffer, BUFFER_SIZE);
                db_seq->clear();
                memset(buffer, 0, BUFFER_SIZE);
                buffer_bit = 0;
            }
        }
        *db_name << name << '\n';
        // insert_seq(db_info, db_i, stream_i, n_ptr, name);
        db_sofs->write((char *)&ptr, sizeof(size_t));
        db_nofs->write((char *)&n_ptr, sizeof(size_t));
        s_num++;
        // *db_sofs << ptr << endl;
        // *db_nofs << n_ptr << endl;
        ptr += content.size();
        n_ptr += name.size() + 1;
    }
    // *db_seq << buffer;
    while (buffer_bit < BUFFER_SIZE * 8)
    {
        buffer_insert(buffer, buffer_bit, -1, content_gaps);
    }
    db_seq->write(buffer, BUFFER_SIZE);
    memset(buffer, 0xFF, BUFFER_SIZE);
    db_seq->write(buffer, BUFFER_SIZE);
    buffer_bit = 0;
    db_sofs->write((char *)&ptr, sizeof(size_t));
    db_nofs->write((char *)&n_ptr, sizeof(size_t));
    // *db_sofs << ptr << endl;
    // *db_nofs << n_ptr << endl;
    db_seq->close();
    db_seq->clear();
    delete db_seq;
    db_name->close();
    db_name->clear();
    delete db_name;
    // sqlite3_exec(db_info, "commit;", 0, 0, 0);
    db_sofs->close();
    db_sofs->clear();
    delete db_sofs;
    db_nofs->close();
    db_nofs->clear();
    delete db_nofs;
    cout << "100%\t Create db " << db_i << "_" << stream_i << " successful! Num. sequence = " << s_num << ", Sequence length = " << ptr << endl;
    memset(buffer, 0, BUFFER_SIZE);
    input.clear();
    cout << "Re-create last db .." << endl;

    size_t last_db_size = 0;
    for (int s = 0; s <= stream_i; s++)
    {
        string fn = f_name + "_" + to_string(db_i) + "_" + to_string(s) + ".seq";
        ifstream in(fn);
        // cout <<fn<<endl;
        in.seekg(0, ios::end);
        int fsz = ((int)in.tellg() / 5) * 8;
        // cout<<fsz<<endl;
        last_db_size += fsz;
        in.close();
        in.clear();
    }

    size_t each_stream_size = (last_db_size - 1) / NUM_STREAM + 1;
    // cout<<last_db_size<<" "<<each_stream_size<<endl;
    input.seekg(db_begin, ios::beg);

    // delete_info_db(db_info, db_i);

    stream_i = 0;
    ptr = 0;
    n_ptr = 0;
    s_num = 0;
    db_seq = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".seq");
    db_name = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".name");
    db_sofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".sofs");
    db_nofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".nofs");
    // sqlite3_exec(db_info, "begin;", 0, 0, 0);
    name.clear();
    content.clear();

    while (std::getline(input, line).good())
    {
        if (line.empty() || line[0] == '>' || line[0] == '#')
        {
            if (!name.empty())
            {
                content.push_back(-1);
                content_gaps = 0;
                for (int j = 0; j < content.size(); j++)
                {
                    char ch = content[j];
                    buffer_insert(buffer, buffer_bit, ch, content_gaps);
                    if (buffer_bit >= BUFFER_SIZE * 8)
                    {
                        assert(buffer_bit == BUFFER_SIZE * 8);
                        // *db_seq << buffer;
                        db_seq->write(buffer, BUFFER_SIZE);
                        db_seq->clear();
                        memset(buffer, 0, BUFFER_SIZE);
                        buffer_bit = 0;
                    }
                }
                *db_name << name << '\n';
                // insert_seq(db_info, db_i, stream_i, n_ptr, name);
                db_sofs->write((char *)&ptr, sizeof(size_t));
                db_nofs->write((char *)&n_ptr, sizeof(size_t));
                s_num++;
                // *db_sofs << ptr << endl;
                // *db_nofs << n_ptr << endl;
                ptr += content.size() - content_gaps;
                n_ptr += name.size() + 1;
                name.clear();
                if (ptr >= each_stream_size)
                {
                    // *db_seq << buffer;
                    while (buffer_bit < BUFFER_SIZE * 8)
                    {
                        buffer_insert(buffer, buffer_bit, -1, content_gaps);
                    }
                    db_seq->write(buffer, BUFFER_SIZE);
                    memset(buffer, 0xFF, BUFFER_SIZE);
                    db_seq->write(buffer, BUFFER_SIZE);
                    buffer_bit = 0;
                    db_sofs->write((char *)&ptr, sizeof(size_t));
                    db_nofs->write((char *)&n_ptr, sizeof(size_t));
                    db_seq->close();
                    db_seq->clear();
                    delete db_seq;
                    db_name->close();
                    db_name->clear();
                    delete db_name;
                    // sqlite3_exec(db_info, "commit;", 0, 0, 0);
                    db_sofs->close();
                    db_sofs->clear();
                    delete db_sofs;
                    db_nofs->close();
                    db_nofs->clear();
                    delete db_nofs;
                    cout << (int)((double)input.tellg() / f_size * 100) << "%\t Re-create db " << db_i << "_" << stream_i << " successful! Num. sequence = " << s_num << ", Sequence length = " << ptr << endl;

                    stream_i++;
                    db_seq = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".seq");
                    db_name = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".name");
                    db_sofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".sofs");
                    db_nofs = new ofstream(f_name + "_" + to_string(db_i) + "_" + to_string(stream_i) + ".nofs");
                    // sqlite3_exec(db_info, "begin;", 0, 0, 0);
                    ptr = 0;
                    n_ptr = 0;
                    s_num = 0;
                    memset(buffer, 0, BUFFER_SIZE);
                }
            }

            if (!line.empty())
            {
                if (line[0] == '>')
                {
                    name = line.substr(1);
                }
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
    {
        content.push_back(-1);
        content_gaps = 0;
        for (int j = 0; j < content.size(); j++)
        {
            char ch = content[j];
            buffer_insert(buffer, buffer_bit, ch, content_gaps);
            if (buffer_bit >= BUFFER_SIZE * 8)
            {
                assert(buffer_bit == BUFFER_SIZE * 8);
                // *db_seq << buffer;
                db_seq->write(buffer, BUFFER_SIZE);
                db_seq->clear();
                memset(buffer, 0, BUFFER_SIZE);
                buffer_bit = 0;
            }
        }
        *db_name << name << '\n';
        // insert_seq(db_info, db_i, stream_i, n_ptr, name);
        db_sofs->write((char *)&ptr, sizeof(size_t));
        db_nofs->write((char *)&n_ptr, sizeof(size_t));
        s_num++;
        // *db_sofs << ptr << endl;
        // *db_nofs << n_ptr << endl;
        ptr += content.size();
        n_ptr += name.size() + 1;
    }
    // *db_seq << buffer;
    while (buffer_bit < BUFFER_SIZE * 8)
    {
        buffer_insert(buffer, buffer_bit, -1, content_gaps);
    }
    db_seq->write(buffer, BUFFER_SIZE);
    memset(buffer, 0xFF, BUFFER_SIZE);
    db_seq->write(buffer, BUFFER_SIZE);

    buffer_bit = 0;
    db_sofs->write((char *)&ptr, sizeof(size_t));
    db_nofs->write((char *)&n_ptr, sizeof(size_t));
    // *db_sofs << ptr << endl;
    // *db_nofs << n_ptr << endl;
    db_seq->close();
    db_seq->clear();
    delete db_seq;
    db_name->close();
    db_name->clear();
    delete db_name;
    // sqlite3_exec(db_info, "commit;", 0, 0, 0);
    db_sofs->close();
    db_sofs->clear();
    delete db_sofs;
    db_nofs->close();
    db_nofs->clear();
    delete db_nofs;
    memset(buffer, 0, BUFFER_SIZE);
    assert(stream_i == NUM_STREAM - 1);
    cout << "100%\t Re-create db " << db_i << "_" << stream_i << " successful! Num. sequence = " << s_num << ", Sequence length = " << ptr << endl;
    input.close();
    input.clear();

    // cout << "Building index .."<< endl;

    // create_index(db_info);

    // close_db(db_info);
}

int main(int argc, const char **argv)
{
    help(argc, argv);
    if (argc <= 3)
    {
        std::cerr << "usage: " << argv[0] << " <protein_db.fasta> <output db> <batch size (GB)>" << std::endl;
        return -1;
    }

    size_t stream_max_size = (size_t)4 * 1024 * 1024 * 1024 - MAX_QUERY_LENGTH;

    size_t db_max_size = stod(argv[3]) * 1024 * 1024 * 1024;

    stream_max_size = min((db_max_size - 1) / NUM_STREAM + 1, stream_max_size);

    create_db(argv[1], argv[2], stream_max_size);
    cout << "Create Database " << argv[2] << " successful! " << endl;
    return 0;
}