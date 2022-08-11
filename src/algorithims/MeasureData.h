#ifndef MEASUREDATA_H
#define MEASUREDATA_H
#include "transform.h"
#include "Point.h"
#include <ctime>
#include <string>
#include <sstream>

using namespace std;



// 测量数据
class MeasureData{
private:
	string ID;						//ID标识符
	double x;					//经度
	double y;					//纬度
	double alpha;				//示向度  360
	double RF;					//RF值
	vector<double> PRI;					//PRI值
	double PW;					//PW值
	struct tm t;				//时间
public:
	//构造函数
	MeasureData()=default;
	MeasureData(double x, double y, double alpha, double RF, vector<double> PRI, double PW, struct tm t);

	//set
	void set(double x, double y, double alpha, double RF, vector<double> PRI, double PW, struct tm t);
	void setID(string ID);
	void setPRI(vector<double> PRI);
	void set_xy(double x, double y, double alpha);

	//get
	double get_x();
	double get_y();
	double get_alpha();
	double get_RF();
	vector<double> get_PRI();
	double get_PW();
	struct tm get_time();
	string get_ID();
	//op
	void angle_hu();
	void hu_angle();
};

//测量数据存储
class MeasureDataSet{
private:
	vector<MeasureData> datas;		//测量数据
	int count;						//数目
	int ID;
	bool isequal(MeasureData m1, MeasureData m2);	// 判断两个数据是否相等
	Point center;
public:
	MeasureDataSet();
	MeasureDataSet(vector<MeasureData> datas);

	//读文件
	void read_file(string path);
	void read_file_only_PRI(string path);
	void read_file_only_xy(string path);
	void read_file_only_RFPWPRI(string path);


	//set
	void set(vector<MeasureData> datas);
	void setID(int ID);
	void set_center(Point center);

	//get
	int getID();
	vector<MeasureData> get_data();
	int get_count();
	Point get_center();

	//op
	int find(MeasureData data);		//如果有该数据 返回id
	bool add(MeasureData data);		//如果表里没有该数据 加上
	bool remove(MeasureData data);		//如果表里有该数据  删除
	void change2angle();			//把原本是弧度的转换成度数
};

string to_str_MeasureData(MeasureData md);
string to_str_MeasureDataSet(MeasureDataSet mds);
string to_str_MeasureDataSet(MeasureDataSet mds, int num);
string to_str_vec_MeasureDataSet(vector<MeasureDataSet> vec_mds);
string to_str_vec_MDCenter(vector<MeasureDataSet> vec_mds);

#endif