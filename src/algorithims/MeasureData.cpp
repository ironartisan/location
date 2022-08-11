#include "MeasureData.h"
using namespace std;
/*
	MeasureData
*/


MeasureData::MeasureData(double x_in, double y_in, double alpha_in, double RF_in, vector<double> PRI_in, double PW_in, struct tm t_in){
	x = x_in;
	y = y_in;
	alpha = alpha_in;
	RF = RF_in;
	PRI = PRI_in;
	PW = PW_in;
	t = t_in;
}

//set
void MeasureData::set(double x_in, double y_in, double alpha_in, double RF_in, vector<double> PRI_in, double PW_in, struct tm t_in){
	x = x_in;
	y = y_in;
	alpha = alpha_in;
	RF = RF_in;
	PRI = PRI_in;
	PW = PW_in;
	t = t_in;
}

void MeasureData::setID(string ID_in){
	ID = ID_in;
}

void MeasureData::setPRI(vector<double> PRI_in){
	PRI = PRI_in;
}

void MeasureData::set_xy(double x_in, double y_in, double alpha_in){
	x = x_in;
	y = y_in;
	alpha = alpha_in;
}

//get
double MeasureData::get_x(){ return x; }

double MeasureData::get_y(){ return y; }

double MeasureData::get_RF(){ return RF; }

vector<double> MeasureData::get_PRI(){ return PRI; }

double MeasureData::get_PW(){ return PW; }

double MeasureData::get_alpha(){ return alpha; }

string MeasureData::get_ID(){ return ID; }

struct tm MeasureData::get_time(){ return t; }

//op
void MeasureData::angle_hu(){
	alpha = angle2hu(alpha);
}

void MeasureData::hu_angle(){
	alpha = hu2angle(alpha);
}

/*
	MeasureDataSet
*/
//辅助函数
bool MeasureDataSet::isequal(MeasureData md1, MeasureData md2){
	double double_e = 1e-5;
	if (abs(md1.get_x() - md2.get_x()) < double_e&&
		abs(md1.get_y() - md2.get_y()) < double_e&&
		abs(md1.get_alpha() - md2.get_alpha()) < double_e&&
		abs(md1.get_RF() - md2.get_RF()) < double_e&&
		abs(md1.get_PW() - md2.get_PW()) < double_e&&
		sub_time(md1.get_time(), md2.get_time()) == 0){
		bool flag = true;
		if (md1.get_PRI().size() != md2.get_PRI().size()){
			return true;
		}
		else{
			for (int i = 0; i < md1.get_PRI().size(); i++){
				if (abs(md1.get_PRI().at(i) - md2.get_PRI().at(i)) < double_e){
					return true;
				}
			}
		}
	}
	return false;
	
}

//构造函数
MeasureDataSet::MeasureDataSet(){
	count = 0;
}

MeasureDataSet::MeasureDataSet(vector<MeasureData> datas_in){
	datas = datas_in;
	count = datas.size();
}

//读文件

void MeasureDataSet::read_file(string path){
	/*
	vector<MeasureData> res;
	ifstream myfile(path);		
	while(!myfile.eof()){
		string ID;
		double x;
		double y;
		double alpha;
		double RF;
		double PW;
		double PRI;
		int year, month, date, hour, min, sec;
		myfile >> ID >> x >> y >> alpha >> RF >> PW >> PRI;
		myfile >> year >> month >> date >> hour >> min >> sec;
		MeasureData tem(x, y, alpha, RF, PRI, PW, create_time(year, month, date, hour, min, sec));
		tem.setID(ID);
		res.push_back(tem);
	}
	res.erase(res.begin() + res.size() - 1);
	set(res);
	*/
}

void MeasureDataSet::read_file_only_PRI(string path){
	vector<MeasureData> res;
	ifstream myfile(path);
	string line;
	if (myfile){
		while (getline(myfile, line)){
			MeasureData md;
			stringstream in(line);
			string id;
			in >> id;
			md.setID(id);
			double x = 0;
			double y = 0;
			double alpha = 0;
			double RF = 0;
			double PW = 0;
			int year = 2021, month = 9, date = 6, hour = 15, min = 46, sec = 10;
			double tem_pri;
			vector<double> pri;
			while (in >> tem_pri){
				pri.push_back(tem_pri);
			}
			md.set(x, y, alpha, RF, pri, PW, create_time(year, month, date, hour, min, sec));
			res.push_back(md);
		}
	}
	datas = res;
}

void MeasureDataSet::read_file_only_xy(string path){
	vector<MeasureData> res;
	ifstream myfile(path);
	while (!myfile.eof()){
		string ID;
		double x;
		double y;
		double alpha;
		myfile >> ID >> x >> y >> alpha ;
		vector<double> pri;
		MeasureData tem(x, y, alpha, 0.0, pri, 0.0, create_time(2021,9,13,19,12,20));
		tem.setID(ID);
		res.push_back(tem);
	}
	res.erase(res.begin() + res.size() - 1);

	set(res);
	change2angle();
}

void MeasureDataSet::read_file_only_RFPWPRI(string path){
	vector<MeasureData> res;
	ifstream myfile(path);
	string line;
	if (myfile){
		while (getline(myfile, line)){
			MeasureData md;
			stringstream in(line);
			string id;
			in >> id;
			md.setID(id);
			double x = 0;
			double y = 0;
			double alpha = 0;
			double RF = 0;
			double PW = 0;
			in >> RF >> PW;
			int year = 2021, month = 9, date = 6, hour = 15, min = 46, sec = 10;
			double tem_pri;
			vector<double> pri;
			while (in >> tem_pri){
				pri.push_back(tem_pri);
			}
			md.set(x, y, alpha, RF, pri, PW, create_time(year, month, date, hour, min, sec));
			res.push_back(md);
		}
	}
	datas = res;
}

//set
void MeasureDataSet::set(vector<MeasureData> datas_in){
	datas = datas_in;
	count = datas.size();
}

void MeasureDataSet::setID(int ID_in){
	ID = ID_in;
}

void MeasureDataSet::set_center(Point center_in){
	center = center_in;
}

//get
vector<MeasureData> MeasureDataSet::get_data(){
	return datas;
}

int MeasureDataSet::getID(){
	return ID;
}

int MeasureDataSet::get_count(){
	return count;
}

Point MeasureDataSet::get_center(){
	return center;
}

//op
int MeasureDataSet::find(MeasureData data){
	for (int i = 0; i < datas.size(); i++){
		if (isequal(data, datas.at(i))){
			return i;
		}
	}
	return -1;
}

bool MeasureDataSet::add(MeasureData data){

	datas.push_back(data);
	count++;
	return true;
	
}

bool MeasureDataSet::remove(MeasureData data){
	int i = find(data);
	if (i != -1){
		datas.erase(datas.begin() + i);
		count--;
		return true;
	}
	else{
		return false;
	}
}

void MeasureDataSet::change2angle(){
	vector<MeasureData> mdlist = datas;
	for (int i = 0; i < mdlist.size(); i++){
		mdlist[i].hu_angle();
	}
	set(mdlist);
}

/*
	辅助
*/
string to_str_PRI(vector<double> PRI){
	string res = "";
	res += "(";
	for (int i = 0; i < PRI.size(); i++){
		res += to_string(PRI.at(i)) + ",";
	}
	res += ")";
	return res;
}

string to_str_MeasureData(MeasureData md){
	string s;
	s = md.get_ID() + " " +
		to_string(md.get_x()) + " " +
		to_string(md.get_y()) + " " +
		to_string(md.get_alpha()) + " " +
		to_string(md.get_RF()) + " " +
		to_string(md.get_PW()) + " " +
		to_str_PRI(md.get_PRI()) + " " +
		to_str_time(md.get_time());
	return s;
}

string to_str_MeasureDataSet(MeasureDataSet mds){
	vector<MeasureData> data = mds.get_data();
	string res = "";
	for (int i = 0; i < data.size(); i++){
		string s = to_str_MeasureData(data.at(i));
		res += s + "\n";
	}
	return res;
}

string to_str_vec_MeasureDataSet(vector<MeasureDataSet> vec_mds){
	string res = "";
	for (int i = 0; i < vec_mds.size(); i++){
		res += "\n\nclass  " + to_string(i) + "\n";
		res += to_str_MeasureDataSet(vec_mds.at(i));
	}
	return res;
}

string to_str_MeasureDataSet(MeasureDataSet mds, int num){
	vector<MeasureData> data = mds.get_data();
	string res = "";
	vector<int> ids = rand_num(num, data.size());
	for (int i = 0; i < ids.size(); i++){
		string s = to_str_MeasureData(data.at(ids.at(i)));
		res += s + "\n";
	}
	return res;
}

string to_str_vec_MDCenter(vector<MeasureDataSet> vec_mds){
	string res = "";
	for (int i = 0; i < vec_mds.size(); i++){
		Point temp = vec_mds.at(i).get_center();
		res += to_string(temp.x) + " " + to_string(temp.y) + "\n";
	}
	return res;
}

