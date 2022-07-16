#ifndef _SMOOTHCU_H_
#define _SMOOTHCU_H_

#include "basic_sketch/basic_sketch.h"

#include <assert.h>

#include <bits/stdc++.h>
#include "Param.h"
#include "Common/hash.h"
#include "basic_smooth.h"

#define inf 2147483647
#define eps 1e-6


typedef long long DATA_TYPE;

template<typename ID_TYPE>
class Stage1_CU {
public: 
	Stage1_CU(){
	}
	~Stage1_CU(){
	}
	void init(uint32_t memory) {
		cell_num = memory * 1024 / sizeof(uint32_t) / 3 / S;
		for (int i = 0; i < 3; ++i) {
			size[i] = cell_num * 32 / counter_bit[i];
		}
		// TowerSketch = new uint32_t** [S];
		TowerSketch = (uint32_t ***)CALLOC(S, sizeof(uint32_t**));
		for (int i = 0; i < S; ++i) {
			// TowerSketch[i] = new uint32_t* [3];
			TowerSketch[i] = (uint32_t **)CALLOC(3, sizeof(uint32_t *));
			for (int j = 0; j < 3; ++j) {
				// TowerSketch[i][j] = new uint32_t [cell_num];
				TowerSketch[i][j] = (uint32_t *)CALLOC(cell_num, sizeof(uint32_t));
				memset(TowerSketch[i][j], 0, cell_num * sizeof(uint32_t));
			}
		}
	}
	void Destroy() {
		for (int i = 0; i < S; ++i) {
			for (int j = 0; j < 3; ++j) {
				FREE(TowerSketch[i][j]);
			}
			FREE(TowerSketch[i]);
		}
		FREE(TowerSketch);
	}
	void add(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end, uint32_t cur_freq){
        uint32_t cur_bit = cur_freq;
        TowerSketch[k][i][j] &= (~(((1 << (end - start)) - 1) << start));
		TowerSketch[k][i][j] |= (cur_bit << start);
		assert(cur_bit == get(k, i, j, start, end));
    }

    uint32_t get(uint32_t k, uint32_t i, uint32_t j, uint32_t start, uint32_t end) {
		return (TowerSketch[k][i][j] & (((1 << (end - start)) - 1) << start)) >> start;
	}
	
   void insert(ID_TYPE id, uint32_t win) {
        // std::cout<<"insert in CU"<<std::endl;
        uint32_t freq[_NumOfArray] = {};
        uint32_t array_cell[_NumOfArray];
        uint32_t array_res[_NumOfArray];
		for (int i = 0; i < _NumOfArray; ++i) {
			uint32_t index = hash_s(id, i) % size[i];
			uint32_t cell = index * counter_bit[i] / 32;
			uint32_t res = index - cell * 32 / counter_bit[i];

            freq[i] =get(win % S, i, cell, res, res + counter_bit[i]);
            if (freq[i] >= max_range[i]) 
                freq[i] = UINT32_MAX;
            array_cell[i] = cell;
            array_res[i] = res;
		}
        uint32_t cur_bit = *std::min_element(freq, freq + _NumOfArray);
        cur_bit ++;

        for(int i = 0; i < _NumOfArray; i++){
            if(freq[i] < cur_bit)
                add(win % S, i, array_cell[i], array_res[i], array_res[i] + counter_bit[i], cur_bit);
        }

	}
	void query(ID_TYPE id, uint32_t win, uint32_t* c) {
		win = win % S;
		for (int k = S - 1; k >= 0; --k, win = (win + S - 1) % S) {
			c[k] = inf;
			for (int i = 0; i < 3; ++i) {
				uint32_t index = hash_s(id, i) % size[i];
				uint32_t cell = index * counter_bit[i] / 32;
				uint32_t res = index - cell * 32 / counter_bit[i];
				uint32_t temp = get(win, i, cell, res, res + counter_bit[i]);
				if (temp <= max_range[i]) 
					c[k] = MIN(c[k], temp);
			}
			assert(c[k] <= window_size);
		}
	}
	void clear(uint32_t win) {
		for (int i = 0; i < 3; ++i) {
			memset(TowerSketch[win % S][i], 0, cell_num * sizeof(uint32_t));
		}
	}
// protected:
public:
	uint32_t*** TowerSketch;
	uint32_t size[3];
	uint32_t cell_num;
	
};



// template<typename ID_TYPE>
class basic_smooth_CU: public basic_sketch {
public: 

    using basic_sketch::operator new;
    using basic_sketch::operator new[];
    using basic_sketch::operator delete;
    using basic_sketch::operator delete[];

	basic_smooth_CU() {}
	basic_smooth_CU(int argc, basic_sketch_string *argv) : win_cnt(0), last_timestamp(0) {//win_cnt, last_timestamp 构造时直接初始化值
		init_matrix();

		int memory = argv[0].to_int();
		double stage1_mem = memory * stage_ratio;
		double stage2_mem = memory * (1 - stage_ratio);

		stage1 = (Stage1_CU<DATA_TYPE> *)CALLOC(1,sizeof(Stage1_CU<DATA_TYPE>));
		stage1_cells = stage1_mem * 1024 / sizeof(uint32_t) / 3 / S;
		stage1->init(stage1_mem);

		stage2 = (Stage2<DATA_TYPE> *)CALLOC(1, sizeof(Stage2<DATA_TYPE>));
		stage2_cells = stage2_mem * 1024 / sizeof(Slot<DATA_TYPE>) / bucket_size;
		stage2->init(stage2_mem);
	}

	basic_smooth_CU(const basic_sketch_string &s){
		init_matrix();

		size_t tmp = 0;
		const char *ss = s.c_str();

		memcpy(&win_cnt, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(&last_timestamp, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(&stage1_cells, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		stage1->cell_num = stage1_cells;

		stage1->TowerSketch = (uint32_t ***)CALLOC(S, sizeof(uint32_t**));
		for (int i = 0; i < S; ++i) {
			stage1->TowerSketch[i] = (uint32_t **)CALLOC(3, sizeof(uint32_t *));
			for (int j = 0; j < 3; ++j) {
				stage1->TowerSketch[i][j] = (uint32_t *)CALLOC(stage1_cells, sizeof(uint32_t));
				memcpy(stage1->TowerSketch[i][j], ss + tmp, stage1_cells * sizeof(uint32_t));
				tmp += stage1_cells * sizeof(uint32_t);
			}
		}

		memcpy(&stage2_cells, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		stage2->size = stage2_cells;

		stage2->slot = (Slot<DATA_TYPE>**)CALLOC(stage2_cells, sizeof(Slot<DATA_TYPE>*));
		for (uint32_t i = 0; i < stage2_cells; ++i) {
			stage2->slot[i] = (Slot<DATA_TYPE>*)CALLOC(bucket_size, sizeof(Slot<DATA_TYPE>));
			memcpy(stage2->slot[i], ss + tmp, bucket_size * sizeof(Slot<DATA_TYPE>));
			tmp += bucket_size * sizeof(Slot<DATA_TYPE>);
		}
		
	}


	basic_sketch_string *to_string(){

		char *s1 = (char*)CALLOC(3 + S * 3 * stage1_cells, sizeof(uint32_t));
		size_t tmp = 0;

		memcpy(s1 + tmp, &win_cnt, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(s1 + tmp, &last_timestamp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(s1 + tmp, &stage1_cells, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		for(int i = 0; i < S; ++i){
			for(int j = 0; j < 3; ++j){
				memcpy(s1 + tmp, &stage1->TowerSketch[i][j], stage1_cells * sizeof(uint32_t));
				tmp += stage1_cells * sizeof(uint32_t);
			}
		}

		basic_sketch_string *bs1 = new basic_sketch_string(s1, tmp);

		char *s2_config = (char *)CALLOC(1, sizeof(uint32_t));
		tmp = 0;

		memcpy(s2_config, &stage2_cells, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		basic_sketch_string *bs2_config = new basic_sketch_string(s2_config, tmp); 
		
		char *s2 = (char *)CALLOC(stage2_cells * bucket_size, sizeof(Slot<DATA_TYPE>));
		tmp = 0;

		for(uint32_t i = 0; i < stage2_cells; ++i){
			memcpy(s2 + tmp, &stage2->slot[i], bucket_size * sizeof(Slot<DATA_TYPE>));
			tmp += bucket_size * sizeof(Slot<DATA_TYPE>);
		}

		basic_sketch_string *bs2 = new basic_sketch_string(s2, tmp);

		*bs1 = *bs1 + *bs2_config + *bs2;

		return bs1;
	}

	~basic_smooth_CU() {
		stage1->Destroy();
		stage2->Destroy();
		FREE(stage1);
		FREE(stage2);
	}

	basic_sketch_reply *insert(const int &argc, const basic_sketch_string *argv) {//将整个数据集整体放入insert中
		win_cnt = 0;
		last_timestamp = 0;

		basic_sketch_reply *insert_out = new basic_sketch_reply;

		for(int i = 0; i < argc; ++i){
			basic_sketch_string itemid = argv[i];
			long long id = itemid.to_long_long();
			uint32_t timestamp = i;

			// insert_out->push_back(id);
			// insert_out->push_back((long long)win_cnt);

			if (last_timestamp + window_size < timestamp) {//transition
				stage2->check(win_cnt);
				// insert_out->push_back((long long)stage2->result.size());//debug

				for(auto it = stage2->result.begin(); it != stage2->result.end(); it++){
					insert_out->push_back(it->first);
					insert_out->push_back((long long)(it->second));
				}
			    transition();
			}

			if (stage2->query(id) >= 0) {
				int id_index = 0;
				stage2->insert(id, win_cnt);
				// insert_out->push_back(id);//debug
				// insert_out->push_back((long long)id_index);//debug
				continue;
			}
			stage1->insert(id, win_cnt);
			
			if (win_cnt < S - 1) 
				continue;
			uint32_t c[S] = {};
			stage1->query(id, win_cnt, c);

			stage2->push(id, win_cnt, c);

			// insert_out->push_back(id);//debug
		
		}

		//transition
		stage2->check(win_cnt);

		for(auto it = stage2->result.begin(); it != stage2->result.end(); it++){
			insert_out->push_back(it->first);
			insert_out->push_back((long long)(it->second));
		}
		transition();

		return insert_out;
	}
	void transition() {
		// std::cout << "Window " << win_cnt << " ends.\n";
		// basic_sketch_reply *result = new basic_sketch_reply;
		// stage2->check(win_cnt);

		// for(auto it = stage2->result.begin(); it != stage2->result.end(); it++){
		// 	result->push_back(it->first);
		// 	result->push_back((long long)(it->second));
		// }
		stage2->vector_swap();

		win_cnt++;
		stage1->clear(win_cnt);
		stage2->clear(win_cnt);
		last_timestamp += window_size;
	}

	// basic_sketch_reply *query() {
	// 	transition();

	// 	basic_sketch_reply *query_out = new basic_sketch_reply;
	// 	query_out->push_back("OK");
	// 	return query_out;
	// }

	
	static basic_sketch_reply *Insert(void *o, const int &argc, const basic_sketch_string *argv){
		return ((basic_smooth_CU *)o)->insert(argc, argv);
	}

	// static basic_sketch_reply *Flush(void *o){
	// 	return ((SmoothSketch*)o)->query();
	// }

	static int command_num(){return 1;};

	static basic_sketch_string command_name (int index){
		basic_sketch_string tmp[] = {"insert"};
		return tmp[index];
	}

	static basic_sketch_func command(int index){
		basic_sketch_func tmp[] = {(basic_smooth_CU::Insert)};
		return tmp[index];
	}

	static basic_sketch_string class_name(){return "basic_smooth_CU";};

	static int command_type(int index){
		int tmp[] = {0};
		return tmp[index];
	}

	static char *type_name(){return "SMOOTH_CU";}
private: 
	uint32_t win_cnt;
	uint32_t last_timestamp;
	Stage1_CU<DATA_TYPE>* stage1;
	Stage2<DATA_TYPE>* stage2;
	uint32_t stage1_cells;
	uint32_t stage2_cells;
};

#endif