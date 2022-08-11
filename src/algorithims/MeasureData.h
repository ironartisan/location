#ifndef MEASUREDATA_H
#define MEASUREDATA_H
#include "transform.h"
#include "Point.h"
#include <ctime>
#include <string>
#include <sstream>

using namespace std;



// ��������
class MeasureData{
private:
	string ID;						//ID��ʶ��
	double x;					//����
	double y;					//γ��
	double alpha;				//ʾ���  360
	double RF;					//RFֵ
	vector<double> PRI;					//PRIֵ
	double PW;					//PWֵ
	struct tm t;				//ʱ��
public:
	//���캯��
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

//�������ݴ洢
class MeasureDataSet{
private:
	vector<MeasureData> datas;		//��������
	int count;						//��Ŀ
	int ID;
	bool isequal(MeasureData m1, MeasureData m2);	// �ж����������Ƿ����
	Point center;
public:
	MeasureDataSet();
	MeasureDataSet(vector<MeasureData> datas);

	//���ļ�
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
	int find(MeasureData data);		//����и����� ����id
	bool add(MeasureData data);		//�������û�и����� ����
	bool remove(MeasureData data);		//��������и�����  ɾ��
	void change2angle();			//��ԭ���ǻ��ȵ�ת���ɶ���
};

string to_str_MeasureData(MeasureData md);
string to_str_MeasureDataSet(MeasureDataSet mds);
string to_str_MeasureDataSet(MeasureDataSet mds, int num);
string to_str_vec_MeasureDataSet(vector<MeasureDataSet> vec_mds);
string to_str_vec_MDCenter(vector<MeasureDataSet> vec_mds);

#endif