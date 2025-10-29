/*
 * @Author: Chao Ning
 * @Date: 2022-04-13 08:44:51
 * @LastEditTime: 2025-01-31 22:49:32
 * @LastEditors: Chao Ning
 * @Description: process genotype files
 * @FilePath: \gmat64\geno\geno.hpp
 */
#define EIGEN_USE_MKL_ALL  // !must be before include Eigen

#include "geno.hpp"
#include "string_utils.hpp"
#include "EigenMatrix_utils.hpp"
#include "iterator_utils.hpp"

/**
 * @brief Construct a new GENO::GENO object
 * !Calculate the number of iids and sids
 * @param geno_file 
 */
GENO::GENO(string geno_file) {
    _geno_file = geno_file;
    // the number of IDs
    string fam_file = _geno_file + ".fam";
    ifstream fin(fam_file);
    if (!fin.is_open()) {
        spdlog::error("Failed to open the fam file: {}", fam_file);
        throw std::runtime_error("File open error");
    }
    _num_id = 0;
    string one_line;
    while (getline(fin, one_line)) { // read line by line different from fin.getline(char* s,size_t n)
        _num_id++;
    }
    fin.close();
    // the number of SNPs
    string bim_file = _geno_file + ".bim";
    fin.open(bim_file);
    if (!fin.is_open()) {
        spdlog::error("Failed to open the bim file: {}", bim_file);
        throw std::runtime_error("File open error");
    }
    _num_snp = 0;
    while (getline(fin, one_line)) { // read line by line different from fin.getline(char* s,size_t n)
        _num_snp++;
    }
    fin.close();
}

/**
 * @brief 
 * ! get the number of iids
 * @return long long 
 */
long long GENO::get_num_iid() {
    return _num_id;
}

/**
 * @brief 
 * ! get the number of SNPs
 * @return long long 
 */
long long GENO::get_num_sid() {
    return _num_snp;
}


/**
 * @brief 
 * !get the iids
 * @return vector<string> 
 */
vector<string> GENO::iid_vec(){
    std::string fam_file = _geno_file + ".fam";
    ifstream fin;
    fin.open(fam_file.c_str());
    if (!fin.is_open()) {
        spdlog::error("Failed to open the fam file: {}", _geno_file + ".fam");
        throw std::runtime_error("File open error");
    }
    string one_line;
    vector<string> id_in_geno_vec;
    vector<string> tmp;
    while (getline(fin, one_line)) {
        process_line(one_line);
        tmp = split_string(one_line);
        if(tmp.size() != 6){
            spdlog::error("Line: {} should have 6 columns!", one_line);
            throw std::runtime_error("Incorrect .fam file format");
        }
        id_in_geno_vec.push_back(tmp[1]);
    }
    fin.close();
    return id_in_geno_vec;
}

/**
 * @brief 
 * ! get the sids
 * @return vector<string> 
 */
vector<string> GENO::sid_vec(){
    string bim_file = _geno_file + ".bim";
    ifstream fin;
    fin.open(bim_file.c_str());
    if(!fin.is_open()){
        spdlog::error("Failed to open the bim file: {}", _geno_file + ".bim");
        throw std::runtime_error("File open error");
    }
    std::string one_line;
    vector<string> tmp; //match pattern *chro snpID base allele1 allele2
    vector<string> snp_id_vec;
    while (getline(fin, one_line)) {
        process_line(one_line);
        tmp = split_string(one_line);
        if(tmp.size() != 6){
            spdlog::error("Line: {} should have 6 columns!", one_line);
            throw std::runtime_error("Incorrect .bim file format");
        }
        snp_id_vec.push_back(tmp[1]);
    }
    return snp_id_vec;
}



/**
 * @brief 
 * ! get the snp annotation information of given number of SNPs from start snp
 * @param start_snp 
 * @param num_snp_read 
 * @return vector<string> 
 */
vector<string> GENO::snp_anno(long long start_snp, long long num_snp_read){
    // check the boundary
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        spdlog::error("The start SNP or end SNP exceeds the boundary, please check!");
        exit(1);
    }
    string bim_file = _geno_file + ".bim";
    ifstream fin;
    fin.open(bim_file.c_str());
    if(!fin.is_open()){
        spdlog::error("Fail to open the bim file: {}", bim_file);
        exit(1);
    }
    string one_line;
    for(long long i = 0; i < start_snp; i++){ // Skip lines before start SNP
        getline(fin, one_line);
    }
    vector<string> snp_anno_vec;
    long long count_line = 0;
    while (getline(fin, one_line)) {
        process_line(one_line);
        snp_anno_vec.push_back(one_line);
        count_line++;
        if(count_line >= num_snp_read){
            break;
        }
    }
    return snp_anno_vec;
}


/**
 * @brief 
 * ! get the snp annotation information with given snp index vector
 * @param snp_index_vec 
 * @return vector<string> 
 */
vector<string> GENO::snp_anno_by_snp_index(vector <long long> snp_index_vec){
    /*
    set<long long> snp_index_set(snp_index_vec.begin(), snp_index_vec.end());
    if(snp_index_set.size() != snp_index_vec.size()){
        cout << "Duplicate IDs exist in the given SNPs" << endl;
        exit(1);
    }
    */

    string bim_file = _geno_file + ".bim";
    ifstream fin;
    fin.open(bim_file.c_str());
    if(!fin.is_open()){
        spdlog::error("Fail to open the bim file:{}", bim_file);
        exit(1);
    }
    string one_line;
    vector<string> snp_anno_all_vec;
    while (getline(fin, one_line)) { //one_line *chro snpID base allele1 allele2
        process_line(one_line);
        snp_anno_all_vec.push_back(one_line);
    }


    vector<string> snp_anno_vec;
    for(vector<long long>::iterator it=snp_index_vec.begin(); it != snp_index_vec.end(); it++) {
        snp_anno_vec.push_back(snp_anno_all_vec[*it]);
    }
    return snp_anno_vec;
}


/**
 * @brief 
 * ! Get the allele frequence and number of observed genotypes from bed file
 * @param freq_arr 
 * @param nobs_geno_arr 
 */
void GENO::bed_allele_freq(string in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr){
    // the number of bytes to store on SNP
    long long num_byte_for_one_snp = _num_id / 4; //One byte can store 4 SNPs
    if ( _num_id % 4 != 0) {
       num_byte_for_one_snp += 1;
    }
    // read
    ifstream fin;
    fin.open(in_file, std::ios::binary); // Read the in binary
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: {}", in_file);
        exit(1);
    }

    double code_val;  // coded value for SNP genotype
    double freq = 0.0;  // SNP frequence
    long long iid = 0, isnp = 0;
    long long num_missing_geno = 0;
    char x03 = '3' - 48; // binary 00000011
    char* bytes_vec = new char[num_byte_for_one_snp]; // store readed btyes for one individual
    VectorXd snp_arr = VectorXd::Zero(_num_id); //SNP array for one individual
    freq_arr.setZero(_num_snp);
    nobs_geno_arr.setZero(_num_snp);
    fin.read(bytes_vec, sizeof(char) * 3); // Read the first three bytes
    while (fin.read(bytes_vec, sizeof(char) * num_byte_for_one_snp)) { // read the SNP genotypes for one locus
        iid = 0;
        num_missing_geno = 0;
        for (long long i = 0; i < num_byte_for_one_snp; i++) {
            for (long long j = 0; j < 4; j++) {
                code_val = (bytes_vec[i] >> (2 * j))& x03; // move two bits from left to right; two bits contain a SNP
                code_val = (code_val * code_val + code_val) / 6.0;
                if (std::fabs(code_val - 2.0 / 6.0) < 0.00001) { // set missing genotype as 2
                    num_missing_geno += 1;
                    code_val = 2;
                }
                snp_arr(iid) = 2 - code_val;  // missing genotype is translated into 0, not affect mean
                iid++; // move to next individual
                if(iid >= _num_id) break;
            }
        }
        nobs_geno_arr(isnp) = _num_id - num_missing_geno;
        if(std::fabs(nobs_geno_arr(isnp)) < 0.01){
            freq = 0;
        }else{
            freq = snp_arr.sum() / (2*_num_id - 2*num_missing_geno);  // delete the miising allele
        }
        freq_arr(isnp) = freq;
        isnp++; // move to next SNP
    }
    fin.close();
    delete [] bytes_vec;
}




/**
 * @brief 
 * ! Get the allele frequence and number of observed genotypes from dgeno file
 * @param freq_arr 
 * @param nobs_geno_arr 
 * @param missing 
 */
void GENO::dgeno_allele_freq(string in_file, VectorXd& freq_arr, VectorXd& nobs_geno_arr, vector<string> missing_in_geno_vec){
    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        spdlog::error("Fail to open the genotype file:{}", in_file);
        exit(1);
    }
    string line;
    vector <string> tmp;
    long long igeno = 0; // ith geno
    long long nobs_geno = 0;
    double freq = 0;
    char *endptr = new char[1000];
    freq_arr.setZero(_num_snp);
    nobs_geno_arr.setZero(_num_snp);
    while(getline(fin, line)){
        process_line(line);
        tmp = split_string(line);
        if(tmp.size() - _num_id != 1){
            spdlog::error("The number of genotyped individuals for locus " + std::to_string(igeno + 1) + " should be " + std::to_string(_num_id));
            exit(1);
        }
        // genotypic value mean and NOBS
        freq = 0.0;
        nobs_geno = 0;
        for(long long i = 1; i < tmp.size(); i++){
            if(!is_nan(tmp[i], missing_in_geno_vec)){
                freq += string_to_double(tmp[i], endptr);
                nobs_geno++;
            }
        }
        if(nobs_geno == 0){
            freq = 0;
        }else{
            freq /= 2 * nobs_geno;
        }
        freq_arr(igeno) = freq;
        nobs_geno_arr(igeno) = nobs_geno;
        igeno++;
        if(igeno > _num_snp){
            spdlog::error("The number of loci exceeds " + std::to_string(_num_snp));
            exit(1);
        }
    }
    delete[] endptr;
    if(igeno < _num_snp){
        spdlog::error("The number of loci less than " + std::to_string(_num_snp));
        exit(1);
    }
}


/**
 * @brief 
 * ! Get the allele frequences (or means/2) and number of observed genotypes from genotypic file
 * @param input_geno_fmt 0, bed file, 1 dgeno file, 2 cgeno file
 * @param freq_arr frequence array for bed and dgeno file, mean array for cgeno file
 * @param nobs_geno_arr 
 * @param missing 
 */
void GENO::allele_freq(int input_geno_fmt, VectorXd& freq_arr, VectorXd& nobs_geno_arr, vector<string> missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::bed_allele_freq(in_file, freq_arr, nobs_geno_arr);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".dgeno";
        GENO::dgeno_allele_freq(in_file, freq_arr, nobs_geno_arr, missing_in_geno_vec);
    }else if(input_geno_fmt == 2){
        string in_file = _geno_file + ".cgeno";
        GENO::dgeno_allele_freq(in_file, freq_arr, nobs_geno_arr, missing_in_geno_vec);
    }
}





/**
 * @brief 
 * ! Get the iid index in the fam file with given iids
 * @param id_in_need_vec A vector of given iids
 * @return vector<long long> 
 */
vector<long long> GENO::find_fam_index(vector<string> id_in_need_vec){
    vector <string> id_in_geno_vec = GENO::iid_vec(); // id in fam
    vector<long long> index_vec;
    long long index;
    vector<string>::iterator is;
    for(vector<string>::iterator it = id_in_need_vec.begin(); it < id_in_need_vec.end(); it++){
        is = find(id_in_geno_vec.begin(), id_in_geno_vec.end(), *it);
        if(is == id_in_geno_vec.end()){
            spdlog::error(*it + " is not in the fam file");
            exit(1);
        }
        index = distance(id_in_geno_vec.begin(), is);
        index_vec.push_back(index);
    }
    return index_vec;
}

vector<long long> GENO::find_fam_index_pro(vector<string> id_in_need_vec){
    set<string> id_in_need_set(id_in_need_vec.begin(), id_in_need_vec.end());
    vector <string> id_in_geno_vec = GENO::iid_vec(); // id in fam
    set<string> id_in_geno_set(id_in_geno_vec.begin(), id_in_geno_vec.end());
    set<string> tmp_set = set_difference_(id_in_need_set, id_in_geno_set);
    if(!tmp_set.empty()){
        spdlog::info("ERROR! These iids are not in the fam file: ", 0, 0, "");
        for(auto tmp:tmp_set){
            spdlog::info(tmp, 0, 0, " ");
        }
        spdlog::info("");
        exit(1);
    }

    std::map<string, long long> id_in_geno_map;
    long long val = 0;
    for(auto tmp:id_in_geno_vec){
        id_in_geno_map[tmp] = val++;
    }
    
    vector<long long> index_vec;
    for(auto tmp:id_in_need_vec){
        index_vec.push_back(id_in_geno_map[tmp]);
    }
    return index_vec;
}

/**
 * @brief 
 * ! Get the sid index in the bim file with given sids
 * @param snp_in_need_vec A vector of given sids
 * @return vector<long long> 
 */
vector<long long> GENO::find_bim_index(vector<string> snp_in_need_vec){
    vector <string> snp_id_in_geno_vec = GENO::sid_vec(); // id in bim
    vector<long long> index_vec;
    long long index;
    vector<string>::iterator is;
    for(vector<string>::iterator it = snp_in_need_vec.begin(); it < snp_in_need_vec.end(); it++){
        is = find(snp_id_in_geno_vec.begin(), snp_id_in_geno_vec.end(), *it);
        if(is == snp_id_in_geno_vec.end()){
            spdlog::error(*it + " is not in the bim file");
            exit(1);
        }
        index = distance(snp_id_in_geno_vec.begin(), is);
        index_vec.push_back(index);
    }
    return index_vec;
}


vector<long long> GENO::find_bim_index_pro(vector<string> snp_in_need_vec){
    set<string> snp_in_need_set(snp_in_need_vec.begin(), snp_in_need_vec.end());
    vector <string> snp_id_in_geno_vec = GENO::sid_vec(); // id in bim
    set<string> snp_id_in_geno_set(snp_id_in_geno_vec.begin(), snp_id_in_geno_vec.end());
    set<string> tmp_set = set_difference_(snp_in_need_set, snp_id_in_geno_set);
    if(!tmp_set.empty()){

        spdlog::info("ERROR! These iids are not in the bim file: ", 0, 0, "");
        for(auto tmp:tmp_set){
            spdlog::info(tmp, 0, 0, " ");
        }
        spdlog::info("");
        exit(1);
        
    }

    std::map<string, long long> id_in_geno_map;
    long long val = 0;
    for(auto tmp:snp_id_in_geno_vec){
        id_in_geno_map[tmp] = val++;
    }
    
    vector<long long> index_vec;
    for(auto tmp:snp_in_need_vec){
        index_vec.push_back(id_in_geno_map[tmp]);
    }
    return index_vec;
}


/**
 * @brief 
 * ! Check the  sids in bim and genotypic (dgeno or cgeno) file, they should be in the same order
 */
void GENO::check_geno_id_order(string geno_file){
    ifstream fin;
    fin.open(geno_file.c_str());
    if(!fin.is_open()){
        spdlog::error("Fail to open the genotypic file: " + geno_file);
        exit(1);
    }

    vector<string> snp_id_in_bim = GENO::sid_vec();
    string one_line;
    string snp_id;
    long long igeno = 0;
    while (getline(fin, one_line)) {
        std::istringstream is(one_line);
        is>>snp_id;
        if(snp_id_in_bim[igeno] != snp_id){
            spdlog::error("The SNP IDs with the same index " + std::to_string(igeno) + 
                    " do not match in genotypic and bim files: " + snp_id + " " + snp_id_in_bim[igeno]);
            exit(1);
        }
        igeno++;
    }
}



void GENO::check_bed_numChar(){
    std::ifstream fin(_geno_file + ".bed", std::ios::binary | std::ios::ate);
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: " + _geno_file + ".bed");
        exit(1);
    }

    std::streampos fileSize = fin.tellg();
    fin.close();

    long long fileSizeExpected = (_num_id + 3) / 4 * _num_snp;
    if(fileSize != fileSizeExpected + 3){
        spdlog::error("The size of " + _geno_file + ".bed" + " doesn't correspond to the number of iids and SNPs");
        exit(1);
    }
}

/**
 * @brief 
 * ! Read genotypic matrix from bed file
 * @param in_file bed file
 * @param snp_mat_part snp matrix
 * @param maf_arr frequence array
 * @param missing_rate_arr missing rate
 * @param start_snp start snp
 * @param num_snp_read the number of readed snp
 * @param index_vec a vector of index to reorder the row (individuals) of snp matrix
 */
void GENO::read_bed(string in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, long long start_snp, long long num_snp_read, 
                vector<long long> index_vec){
    // check the boundary
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        spdlog::error("The start SNP or end SNP exceeds the boundary, please check");
        exit(1);
    }
    
    // the number of bytes to store on SNP
    long long num_byte_for_one_snp = _num_id / 4; //One byte can store 4 SNPs
    if ( _num_id % 4 != 0) {
       num_byte_for_one_snp += 1;
    }
    // read
    ifstream fin;
    fin.open(in_file, std::ios::binary); // Read the in binary
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }

    char* bytes_vec = new char[num_byte_for_one_snp * num_snp_read]; // store readed btyes for all individuals
    fin.read(bytes_vec, sizeof(char) * 3); // Read the first three bytes
    fin.seekg(start_snp*num_byte_for_one_snp, std::ios::cur); // move to the position for the start SNP
    fin.read(bytes_vec, sizeof(char) * num_byte_for_one_snp * num_snp_read);

    double code_val;  // coded value for SNP genotype
    double freq = 0.0;  // SNP frequence
    char x03 = '3' - 48; // binary 00000011
    VectorXd snp_arr = VectorXd::Zero(_num_id); //SNP array for one individual
    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    snp_mat_part.setZero(_num_id, num_snp_read);

    for(long long isnp = 0; isnp < num_snp_read; isnp++){
        long long iid = 0;
        vector<long long> missing_geno_index;
        for (long long i = 0; i < num_byte_for_one_snp; i++) {
            char one_byte = bytes_vec[i + isnp*num_byte_for_one_snp];
            for (long long j = 0; j < 4; j++) {
                code_val = (one_byte >> (2 * j))& x03; // move two bits from left to right; two bits contain a SNP
                code_val = (code_val * code_val + code_val) / 6.0;
                if (std::fabs(code_val - 2.0 / 6.0) < 0.00001) { // set missing genotype as 2
                    missing_geno_index.push_back(iid);
                    code_val = 2;
                }
                snp_arr(iid) = 2 - code_val;  // missing genotype is translated into 0, not affect mean
                iid++; // move to next individual
                if(iid >= _num_id) break;
            }
        }
        // missing rate
        long long num_missing_geno = missing_geno_index.size();
        missing_rate_arr(isnp) = num_missing_geno * 1.0 / _num_id;
        // maf
        if(std::fabs(missing_rate_arr(isnp)) >= 0.99999){
            freq = 0;
        }else{
            freq = snp_arr.sum() / (2 * _num_id - 2 * num_missing_geno);  // delete the missing genotypes
        }
        // set missing genotypes as mean values
        for(vector<long long>::iterator it = missing_geno_index.begin(); it != missing_geno_index.end(); it++){ 
            snp_arr(*it) = 2 * freq; // mean = 2 *freq
        }
        maf_arr(isnp) = freq;
        snp_mat_part.col(isnp) = snp_arr;
    }
    delete [] bytes_vec;
    // reorder the snp matrix
    if(!index_vec.empty()){
        snp_mat_part = (snp_mat_part(index_vec, Eigen::all)).eval();
        maf_arr = snp_mat_part.colwise().sum() / (2 * index_vec.size());
    }
    fin.close();
}


void GENO::read_bed_omp_spMVLMM(string in_file, double* snp_mat_part_pt, VectorXd& maf_arr, VectorXd& missing_rate_arr, long long start_snp, long long num_snp_read, 
                vector<long long> index_vec){
    // check the boundary
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        spdlog::error("The start SNP or end SNP exceeds the boundary, please check");
        exit(1);
    }
    // Compute the number of bytes needed to store genotypes for all individual per SNP.
    long long num_byte_for_one_snp = (_num_id + 3) / 4; //One byte can store 4 SNPs
    
    // read
    FILE* fin = fopen(in_file.c_str(), "rb"); // Read the in binary
    if (!fin) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }

    // LOGGER.i("Read the PLINK BED file");
    // LOGGER.i("Seek");
    std::streamoff startPos = start_snp*num_byte_for_one_snp + 3; // Calculate the starting position of the SNP data
    if (fseek(fin, startPos, SEEK_SET) != 0){
        spdlog::error("Fail to seek the predefined position");
        exit(1);
    }
    
    // LOGGER.i("Read");
    char* bytes_vec = new char[num_byte_for_one_snp * num_snp_read]; // store readed btyes for all individuals
    long long read_count = fread(bytes_vec, sizeof(char), num_byte_for_one_snp * num_snp_read, fin);
    if (read_count != num_byte_for_one_snp * num_snp_read){
        delete [] bytes_vec;
        spdlog::error("fread failed to read the expected amount");
        exit(1);
    }
    if(ferror(fin)){
        delete [] bytes_vec;
        spdlog::error("Failed to read from file: " + in_file);
        exit(1);
    }
    fclose(fin);

    // LOGGER.i("Decode the genotype");
    char x03 = '3' - 48; // binary 00000011
    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);

    long long num_used_id = _num_id;
    if(!index_vec.empty()){
        num_used_id = index_vec.size();
    }
    long long *index_vec_pt = new long long[num_used_id];

    if(index_vec.empty()){
        for(long long i = 0; i < num_used_id; i++)
            index_vec_pt[i] = i;
    }else{
        long long k = 0;
        for(auto tmp:index_vec) 
            index_vec_pt[k++] = tmp;
    }

    long long _omp_max_threads = omp_get_max_threads();
    double **snp_arr_part = new double*[_omp_max_threads];
    for(int i = 0; i < _omp_max_threads; i++){
        snp_arr_part[i] = new double[_num_id];
    }

    #pragma omp parallel for schedule(dynamic)
    for(long long isnp = 0; isnp < num_snp_read; isnp++){
        
        long long thread_id = omp_get_thread_num();
        double freq = 0; // SNP frequence
        long long num_missing = 0;
        double code_val;

        for(long long i = 0; i < _num_id; i++){
            unsigned char byte = bytes_vec[i / 4 + isnp*num_byte_for_one_snp];
            unsigned char genotype = (byte >> (2 * (i % 4))) & 3;
            
            switch (genotype) {
                case 0: code_val = 2.0; freq += code_val; break;  // Homozygous for the major allele.
                case 1: code_val = -9; num_missing++; break;  // Missing genotype.
                case 2: code_val = 1.0; freq += code_val; break;  // Heterozygous.
                case 3: code_val = 0.0; break;  // Homozygous for the minor allele.
                default: throw std::logic_error("Invalid genotype value");
            }

            snp_arr_part[thread_id][i] = code_val;
            
        }

        // missing rate
        missing_rate_arr(isnp) = num_missing * 1.0 / _num_id;

        // maf
        if(num_missing != _num_id){
            freq /= (2 * _num_id - 2 * num_missing);
        }else{
            freq = 0.0;
        }
        maf_arr(isnp) = freq;
        
        // -2*p, nan->0
        long long code_int;
        for(long long i = 0; i < num_used_id; i++){
            code_int = snp_arr_part[thread_id][index_vec_pt[i]];
            if(code_int == -9){
                code_val = 0;
            }else{
                code_val = code_int - 2 * freq;
            }
            snp_mat_part_pt[isnp*num_used_id + i] = code_val;
        }
        
    }
    
    // free memory
    delete [] bytes_vec;
    delete [] index_vec_pt;
    for(int i = 0; i < _omp_max_threads; i++){
        delete [] snp_arr_part[i];
    }
    delete [] snp_arr_part;
    
    // LOGGER.i("Finish");
}


void GENO::read_bed_omp(string in_file, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, long long start_snp, long long num_snp_read, 
                vector<long long> index_vec){
    // check the boundary
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        spdlog::error("The start SNP or end SNP exceeds the boundary, please check");
        exit(1);
    }
    // the number of bytes to store on SNP
    long long num_byte_for_one_snp = _num_id / 4; //One byte can store 4 SNPs
    if ( _num_id % 4 != 0) {
       num_byte_for_one_snp += 1;
    }
    // read
    ifstream fin;
    fin.open(in_file, std::ios::binary); // Read the in binary
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }

    // spdlog::info("Read the PLINK BED file");
    char* bytes_vec = new char[num_byte_for_one_snp * num_snp_read]; // store readed btyes for all individuals
    fin.read(bytes_vec, sizeof(char) * 3); // Read the first three bytes
    fin.seekg(start_snp*num_byte_for_one_snp, std::ios::cur); // move to the position for the start SNP
    fin.read(bytes_vec, sizeof(char) * num_byte_for_one_snp * num_snp_read);


    // spdlog::info("Decode the genotype");
    char x03 = '3' - 48; // binary 00000011
    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);

    long long num_used_id;
    if(index_vec.empty()){
        num_used_id = _num_id;
    }else{
        num_used_id = index_vec.size();
    }


    double *snp_mat_part_pt = new double[num_snp_read * num_used_id];
    long long *index_vec_pt = new long long[num_used_id];
    long long k = 0;
    for(auto tmp:index_vec) index_vec_pt[k++] = tmp;

    #pragma omp parallel for schedule(dynamic)
    for(long long isnp = 0; isnp < num_snp_read; isnp++){
        long long iid = 0;
        vector<long long> missing_geno_index;
        // VectorXd snp_arr = VectorXd::Zero(_num_id); //SNP array for one individual
        double *snp_arr = new double[_num_id];
        double freq = 0; // SNP frequence
        for (long long i = 0; i < num_byte_for_one_snp; i++) {
            char one_byte = bytes_vec[i + isnp*num_byte_for_one_snp];
            for (long long j = 0; j < 4; j++) {
                double code_val = (one_byte >> (2 * j))& x03; // move two bits from left to right; two bits contain a SNP
                code_val = (code_val * code_val + code_val) / 6.0;
                if (std::fabs(code_val - 2.0 / 6.0) < 0.00001) { // set missing genotype as 2
                    missing_geno_index.push_back(iid);
                    code_val = 2;
                }
                code_val = 2 - code_val;
                snp_arr[iid] = code_val;  // missing genotype is translated into 0, not affect mean
                freq += code_val;
                iid++; // move to next individual
                if(iid >= _num_id) break;
            }
        }
        
        // missing rate
        long long num_missing_geno = missing_geno_index.size();
        missing_rate_arr(isnp) = num_missing_geno * 1.0 / _num_id;

        // maf
        if(std::fabs(missing_rate_arr(isnp)) >= 0.99999){
            freq = 0;
        }else{
            freq /= (2 * _num_id - 2 * num_missing_geno);  // delete the missing genotypes
        }
        maf_arr(isnp) = freq;
        // set missing genotypes as mean values
        for(auto tmp:missing_geno_index){
            snp_arr[tmp] = 2 * freq;
        }


        if(index_vec.empty()){
            for(long long i = 0; i < _num_id; i++){
                snp_mat_part_pt[isnp*_num_id + i] = snp_arr[i];
            }
        }else{
            for(long long i = 0; i < num_used_id; i++){
                snp_mat_part_pt[isnp*num_used_id + i] = snp_arr[index_vec_pt[i]];
            }
        }
        delete [] snp_arr;
    }

    // spdlog::info("Map");
    snp_mat_part = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic > > (snp_mat_part_pt, num_used_id, num_snp_read);
    /*
    snp_mat_part.resize(num_used_id, num_snp_read);
    #pragma omp parallel for schedule(dynamic)
    for(long long i = 0; i < num_used_id; i++){
        for(long long j = 0; j < num_snp_read; j++)
            snp_mat_part(i, j) = snp_mat_part_pt[j*num_used_id + i];
    }
    */
    delete [] bytes_vec;
    delete [] index_vec_pt;
    delete [] snp_mat_part_pt;
    // spdlog::info("Finish");
    fin.close();
}



/**
 * @brief 
 * ! Read genotypic matrix from dgeno file
 * @param in_file dgeno file
 * @param geno_mat_part genotypic matrix
 * @param freq_arr frequence array = mean/2
 * @param missing_rate_arr missing rate
 * @param start_snp start snp
 * @param num_snp_read the number of readed snp
 * @param index_vec a vector of index to reorder the row (individuals) of snp matrix
 * @param missing string used to define missing values
 */
void GENO::read_dgeno(string in_file, MatrixXd& geno_mat_part, VectorXd& freq_arr, VectorXd& missing_rate_arr,
            long long start_snp, long long num_snp_read, vector<long long> index_vec, vector<string> missing_in_geno_vec){
                
    // check the boundary
    if(start_snp < 0 || start_snp + num_snp_read > _num_snp || num_snp_read <= 0){
        spdlog::error("The start SNP or end SNP exceeds the boundary, please check");
        std::exit(1);
    }

    
    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        spdlog::error("Fail to open the genotypic file: " + in_file);
        std::exit(1);
    }

    string line;
    vector <string> tmp;
    // Skip lines before start SNP
    for(long long i = 0; i < start_snp; i++){
        getline(fin, line);
    }

    
    long long num_missing_geno;
    double freq = 0;
    freq_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    geno_mat_part.setZero(_num_id, num_snp_read);
    // Read
    long long igeno = 0;
    VectorXd geno_arr = VectorXd::Zero(_num_id); //SNP array for one individual
    char *endptr = new char[1000];
    while(getline(fin, line)){
        process_line(line);
        tmp = split_string(line);
        if(tmp.size() - _num_id != 1){
            spdlog::error("The number of genotyped individuals for locus " + std::to_string(igeno + 1) + " should be " + std::to_string(_num_id));
            exit(1);
        }
        // genotypic values
        freq = 0.0;
        num_missing_geno = 0;
        vector<long long> missing_geno_index; // the index for missing genotypes
        for(long long i = 1; i < tmp.size(); i++){
            if(!is_nan(tmp[i], missing_in_geno_vec)){
                freq += string_to_double(tmp[i], endptr);
                geno_arr(i - 1) = string_to_double(tmp[i], endptr);
            }else{
                num_missing_geno++;
                missing_geno_index.push_back(i - 1);
            }
        }
        // frequence values and set missing genotypes as mean
        if(num_missing_geno == _num_id){
            freq = 0;
        }else{
            freq /= 2*_num_id - 2*num_missing_geno;
        }
        for(vector<long long>::iterator it = missing_geno_index.begin(); it != missing_geno_index.end(); it++){ // set missing genotypes as mean values
            geno_arr(*it) = freq * 2;
        }
        geno_mat_part.col(igeno) = geno_arr;
        freq_arr(igeno) = freq;
        igeno++; // move to next SNP
        if(igeno >= num_snp_read)
            break;
    }
    delete[] endptr;
    // reorder the snp matrix
    if(!index_vec.empty()){
        geno_mat_part = geno_mat_part(index_vec, Eigen::all).eval();
        freq_arr = geno_mat_part.colwise().sum() / (2 * index_vec.size());
    }
    fin.close();
}


/**
 * @brief 
 * ! read genotypic matrix using different input formats
 * @param input_geno_fmt input format: 0, bed file; 1, dgeno file; 2, cgeno file
 * @param snp_mat_part 
 * @param maf_arr 
 * @param missing_rate_arr 
 * @param start_snp 
 * @param num_snp_read 
 * @param index_vec 
 * @param missing 
 */
void GENO::read_geno(int input_geno_fmt, MatrixXd& snp_mat_part, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
            long long start_snp, long long num_snp_read, vector<long long> index_vec, vector<string> missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::read_bed_omp(in_file, snp_mat_part, maf_arr, missing_rate_arr, start_snp, num_snp_read, index_vec);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".fmat"; // feture matrix such as SNP, gene expression, etc
        GENO::read_dgeno(in_file, snp_mat_part, maf_arr, missing_rate_arr, start_snp, num_snp_read, index_vec, missing_in_geno_vec);
    }else{
        spdlog::error("The input genotypic format is not supported");
        std::exit(1);
    }
}


/**
 * @brief 
 * ! Read genotypic matrix with given snp index
 * @param in_file 
 * @param snp_mat_by_snp_index 
 * @param maf_arr 
 * @param missing_rate_arr 
 * @param snp_index_vec 
 * @param index_vec 
 */
void GENO::read_bed_by_snp_index(string in_file, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                vector<long long> snp_index_vec, vector<long long> index_vec){
    
    // the number of bytes to store on SNP
    long long num_byte_for_one_snp = _num_id / 4; //One byte can store 4 SNPs
    if ( _num_id % 4 != 0) {
       num_byte_for_one_snp += 1;
    }

    // read
    ifstream fin;
    fin.open(in_file, std::ios::binary); // Read the in binary
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }


    double code_val, missing_rate;  // coded value, missing rate
    double freq = 0.0;  // SNP frequence
    long long iid = 0, isnp = 0;
    long long num_missing_geno = 0;
    vector<long long> missing_geno_index; // the index for missing genotypes
    
    char x03 = '3' - 48; // binary 00000011
    char* bytes_vec = new char[num_byte_for_one_snp]; // store readed btyes for one individual
    VectorXd snp_arr = VectorXd::Zero(_num_id); //SNP array for one individual
    long long num_snp_read = snp_index_vec.size();
    maf_arr.setZero(num_snp_read);
    missing_rate_arr.setZero(num_snp_read);
    snp_mat_by_snp_index.setZero(_num_id, num_snp_read);
    fin.read(bytes_vec, sizeof(char) * 3); // Read the first three bytes
    for(vector<long long>::iterator it=snp_index_vec.begin(); it != snp_index_vec.end(); it++) { // read the SNP genotypes for one locus
        fin.clear();
        iid = 0;
        num_missing_geno = 0;
        missing_geno_index.clear();
        isnp = distance(snp_index_vec.begin(), it); // snp index
        // seek the position to read
        fin.seekg(*it*num_byte_for_one_snp + 3, std::ios::beg);
        fin.read(bytes_vec, sizeof(char) * num_byte_for_one_snp);
        for (long long i = 0; i < num_byte_for_one_snp; i++) {
            for (long long j = 0; j < 4; j++) {
                code_val = (bytes_vec[i] >> (2 * j))& x03; // move two bits from left to right; two bits contain a SNP
                code_val = (code_val * code_val + code_val) / 6.0;
                if (std::fabs(code_val - 2.0 / 6.0) < 0.00001) { // set missing genotype as 2
                    num_missing_geno++;
                    missing_geno_index.push_back(iid);
                    code_val = 2;
                }
                snp_arr(iid) = 2 - code_val;  // missing genotype is translated into 0, not affect mean
                iid++; // move to next individual
                if(iid >= _num_id) break;
            }
        }
        // missing rate, allele mean
        missing_rate = num_missing_geno / _num_id;
        missing_rate_arr(isnp) = missing_rate;
        if(missing_rate >= 0.99999){
            freq = 0;
        }else{
            freq = snp_arr.sum() / (2 * _num_id - 2 * num_missing_geno);  // delete the missing allele
        }
        // set missing values as allele mean
        for(vector<long long>::iterator it2 = missing_geno_index.begin(); it2 != missing_geno_index.end(); it2++){ // set missing genotypes as mean values
            snp_arr(*it2) = 2 * freq; // mean = 2 * p
        }
        maf_arr(isnp) = freq;
        snp_mat_by_snp_index.col(isnp) = snp_arr;
    }
    fin.close();
    delete [] bytes_vec;
    if(!index_vec.empty()){
        snp_mat_by_snp_index = (snp_mat_by_snp_index(index_vec, Eigen::all)).eval();
        maf_arr = snp_mat_by_snp_index.colwise().sum() / (2 * index_vec.size());
    }
}


/**
 * @brief 
 * ! Read genotypic matrix with given snp index
 * @param in_file 
 * @param snp_mat_by_snp_index 
 * @param maf_arr 
 * @param missing_rate_arr 
 * @param snp_index_vec 
 * @param index_vec 
 * @param missing missing values
 */
void GENO::read_dgeno_by_geno_index(string in_file, MatrixXd& geno_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                vector<long long> geno_index_vec, vector<long long> index_vec, vector<string> missing_in_geno_vec){
    
    // open
    ifstream fin;
    fin.open(in_file);
    if (!fin.is_open()) {
        spdlog::error("Fail to open the dosage genotype file: " + in_file);
        exit(1);
    }


    string line;
    vector <string> tmp;
    long long num_missing_geno;
    double freq = 0;
    long long num_geno_read = geno_index_vec.size();
    maf_arr.setZero(num_geno_read);
    missing_rate_arr.setZero(num_geno_read);
    geno_mat_by_snp_index.setZero(_num_id, num_geno_read);

    // Read
    long long igeno = 0;
    char *endptr = new char[1000];
    VectorXd geno_arr = VectorXd::Zero(_num_id); //genotype array for one individual
    vector<long long>::iterator it;
    long long finded_position;
    while(getline(fin, line)){

        it = find(geno_index_vec.begin(), geno_index_vec.end(), igeno); // locate the needed line
        if(it == geno_index_vec.end()){
            igeno++;
            continue;
        }else{
            finded_position = distance(geno_index_vec.begin(), it);
            // finded_position = &*it-&snp_index_vec[0];
        }

        process_line(line);
        tmp = split_string(line);
        if(tmp.size() - _num_id != 1){
            spdlog::error("The number of genotyped individuals for locus " + std::to_string(igeno + 1) + " should be " + std::to_string(_num_id));
            exit(1);
        }

        
        // genotypic values
        freq = 0.0;
        num_missing_geno = 0;
        vector<long long> missing_geno_index; // the index for missing genotypes
        for(long long i = 1; i < tmp.size(); i++){
            if(!is_nan(tmp[i], missing_in_geno_vec)){
                freq += string_to_double(tmp[i], endptr);
                geno_arr(i - 1) = string_to_double(tmp[i], endptr);
            }else{
                num_missing_geno++;
                missing_geno_index.push_back(i - 1);
            }
        }

        // freq
        if(num_missing_geno == _num_id){
            freq  = 0;
        }else{
            freq /= 2 * _num_id - 2 * num_missing_geno;
        }
        maf_arr(finded_position) = freq;


        // set missing genotypes as mean values
        for(vector<long long>::iterator it = missing_geno_index.begin(); it != missing_geno_index.end(); it++){ 
            geno_arr(*it) = 2 * freq;
        }
        
        //Snp matrix
        geno_mat_by_snp_index.col(finded_position) = geno_arr;
        igeno++; // move to next SNP
    }
    delete[] endptr;
    fin.close();
    if(!index_vec.empty()){
        geno_mat_by_snp_index = geno_mat_by_snp_index(index_vec, Eigen::all).eval();
        maf_arr = geno_mat_by_snp_index.colwise().sum() / (2 * index_vec.size());
    }
}


/**
 * @brief 
 * ! Read genotypic matrix with given snp index from different input genotypic formats
 * @param input_geno_fmt 
 * @param snp_mat_by_snp_index 
 * @param maf_arr 
 * @param missing_rate_arr 
 * @param snp_index_vec 
 * @param index_vec 
 * @param missing 
 */
void GENO::read_geno_by_index(int input_geno_fmt, MatrixXd& snp_mat_by_snp_index, VectorXd& maf_arr, VectorXd& missing_rate_arr, 
                vector<long long> snp_index_vec, vector<long long> index_vec, vector<string> missing_in_geno_vec){
    if(input_geno_fmt == 0){
        string in_file = _geno_file + ".bed";
        GENO::read_bed_by_snp_index(in_file, snp_mat_by_snp_index, maf_arr, missing_rate_arr, 
                snp_index_vec, index_vec);
    }else if(input_geno_fmt == 1){
        string in_file = _geno_file + ".dgeno";
        GENO::read_dgeno_by_geno_index(in_file, snp_mat_by_snp_index, maf_arr, missing_rate_arr, 
                snp_index_vec, index_vec, missing_in_geno_vec);
    }else{
        string in_file = _geno_file + ".cgeno";
        GENO::read_dgeno_by_geno_index(in_file, snp_mat_by_snp_index, maf_arr, missing_rate_arr, 
                snp_index_vec, index_vec, missing_in_geno_vec);
    }
}


/**
 * @brief
 * ! Read bed file with start SNP index and iid index
 * @param snp_mat_part snp matrix
 * @param start_snp 
 * @param num_snp_read 
 * @param start_id 
 * @param num_id_read 
 */
void GENO::read_bed_2Dpart(MatrixXd& snp_mat_part, long long start_snp, long long num_snp_read, 
        long long start_id, long long num_id_read){
    // the number of bytes to store on SNP
    long long num_byte_for_one_snp = _num_id / 4; //One byte can store 4 SNPs
    if ( _num_id % 4 != 0) {
       num_byte_for_one_snp += 1;
    }
    string in_file = _geno_file + ".bed";
    ifstream fin;
    fin.open(in_file, std::ios::binary); // Read the in binary
    if (!fin.is_open()) {
        spdlog::error("Fail to open the plink bed file: " + in_file);
        exit(1);
    }
    char* bytes_mat = new char[num_byte_for_one_snp * num_snp_read]; // store readed btyes for all individuals
    double code_val, mean_allele;  // coded value for SNP genotype
    long long iid = 0, isnp = 0;
    long long num_missing_allele = 0;
    char x03 = '3' - 48; // binary 00000011
    VectorXd snp_arr = VectorXd::Zero(_num_id); //SNP array for one individual, ordered same as fam file
    fin.read(bytes_mat, sizeof(char) * 3); // Read the first three bytes
    fin.seekg(start_snp*num_byte_for_one_snp, std::ios::cur); // move to the position for the start SNP
    fin.read(bytes_mat, sizeof(char) * num_byte_for_one_snp * num_snp_read);
    snp_mat_part.setZero(num_id_read, num_snp_read);
    for(long long k = 0; k < num_snp_read; k++) { // read the SNP genotypes for all individuals
        vector<long long> missing_geno_index; // the index for missing genotypes
        for (long long i = 0; i < num_byte_for_one_snp; i++) {
            for (long long j = 0; j < 4; j++) {
                code_val = (bytes_mat[k * num_byte_for_one_snp + i] >> (2 * j))& x03; // move two bits from left to right; two bits contain a SNP
                code_val = (code_val * code_val + code_val) / 6.0;
                if (std::fabs(code_val - 2.0 / 6.0) < 0.00001) { // missing genotype set as 2
                    num_missing_allele += 2;
                    missing_geno_index.push_back(iid);
                    code_val = 2;
                }
                snp_arr(iid++) = 2 - code_val;  // missing genotype is translated into 0, not affect mean
                if(iid >= _num_id) break;
            }
        }
        if(num_missing_allele == 2 * _num_id){
            mean_allele = 0;
        }else{
            mean_allele = snp_arr.sum() / (_num_id - num_missing_allele/2);  // delete the missing allele
        }
        // set missing genotypes as mean values
        for(vector<long long>::iterator it = missing_geno_index.begin(); it != missing_geno_index.end(); it++){ 
            snp_arr(*it) = mean_allele;
        }
        snp_mat_part.col(k) = snp_arr.segment(start_id, num_id_read);
        iid = 0;
    }
    delete [] bytes_mat;
    fin.close();
}
