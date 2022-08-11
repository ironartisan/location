#ifndef PLandEE_H
#define PLandEE_H
#include<vector>
#include<cstdio>
#include<iostream>
#include<cmath>
#include<ctime>
#include<algorithm>
#include "Point.h"
#include "Kmeans.h"

using namespace std;
//����ת��
vector<double> GeographicToGeodeticRectangular(double Long, double Lati, double Height);
vector<double> GeodeticRectangularToStationRectangular(double xg, double yg, double zg, double Long, double Lati, double Height);

vector<double> StationRectangularToGeodeticRectangular(double xStationTarget, double yStationTarget, double zStationTarget, double longitudeStation, double latitudeStation, double heightStation);
vector<double> GeodeticRectangularToGeographic(double dX, double dY, double dZ);

double PicAngle2AlgAngle(double in);
double AlgAngle2PicAngle(double in);

// DF data �������� Ҳ����ԭʼ������
class DFdata{
private:
	double x; //����
	double y; //γ��
	double z; //�߶�
	double alpha; //ƽ���ƫ��  0  --- 2pi
	double beta;  //�߶�ƫ��   -pi --- pi
	double r[2];	  // ��ǰ����������ԭ��ľ���
	double sk[2];	  // �����ĴӲ����㵽��λ��ĵ�λ����
	double distance; // �붨λ��ľ��� xyƽ���ϵ�
public:
	DFdata();
	DFdata(double x_in, double y_in, double z_in, double alpha_in, double beta_in, bool isPI); // in  input����˼  

	//���� ��ȡ�͸���
	void set_x(double x_in);
	void set_y(double y_in);
	void set_z(double z_in);
	double set_distance(double x_in, double y_in);	//����Ķ�λ�� �� ������ľ���
	double get_distance();
	// ���ڽǶ���������ã� �Ƕ���������λ ���� �� �Ƕ� �����ʱ��˵���ǻ��Ȼ��ǽǶ� isPI�ǻ���  isPI = 1 ����  isPI = 0 �Ƕ�
	void set_alpha(double alpha_in, bool isPI);
	void set_beta(double beta_in, bool isPI);
	double get_alpha();
	double get_beta();
	vector<double> get_point();  //3*1
	// �й��㷨�в���r��ʽ��  r 2*1 ��ά����
	vector<double> set_r(double x_in, double y_in);
	vector<double> get_r();
	// ����sk��غ��� sk 2*1 ��ά���� skԭ���ǴӲ����� �� ��λ�㴹ֱ�ڲ������ϵĵ㣬 ������Ϊ��λ�㲻�̶������������ò����߷���ĵ�λ����
	vector<double> get_sk();
	vector<double> get_skt();
	// �ж��Ƿ���  ���� ��ֱ��
	bool isValid(double x, double y);
};

// ���ڲ������ݵ���Դ��λ�ⷨ  Q1 �㷨һ
class PLocation{
private:
	vector<DFdata> dflist; // ͬһ���ε�ԭʼ����
	vector<vector<double>> A;		// ����A skt  N*2 
	vector<double> b;		// ����b   skt * r  N*1
	int n;					//��һ���ε����ݵ�����
	double x;					//��λ�㾭��
	double y;					//��λ��ά��
	double z;					//��λ��߶�
	double Pe;					//�����Բ׼ȷ��
	double rl;					//���� r long
	double rs;					//���� r short
	double fan;					//ƫ��
	double sig11;
	double sig12;
	double sig22;
	bool Valid;					//�Ƿ���Ч
	vector<double> stationloc;	//����վ��
	vector<vector<double>> result;		//������
	//������ǰ�㷨����ֵ
	vector<vector<double>> pinvA; //pinv(A) ���� (ATA)^-1 * A    2 * N
	double help_b(vector<double> skt, vector<double> r); //��������b�Ĳ���
	// �㷨��ִ�в���
	vector<vector<double>> get_ATA(); //���AT * A  2*2 �ľ���
	bool ATAsuit(vector<vector<double>> ATA); //�ж�ATA�Ƿ�Ϸ�
	vector<vector<double>> get_reATA(vector<vector<double>> ATA);	//	��ȡATA����
	vector<vector<double>> get_pinvA();//��ȡ�㷨�в���pinv(A)  �� (ATA)^-1 * A    2 * N�ľ���
	vector<double> set_locxy(); // ��ȡxyλ��  
	void set_ellipse();			// �������Բ	
	vector<vector<double>> get_A();	//��ȡA����
	vector<double> get_b(); //��ȡb����
	double set_z();				// ���ݶ�λ����z�Ĺ���ֵ 
	double func_fan(double sig11, double sig12, double sig22);	//����x�� y�� ����õ�ƫ��.


public:
	PLocation();			//�޲ι��캯��
	PLocation(vector<DFdata> dflist_in);	//�вι��캯��
	void input(vector<double> x_in, vector<double> y_in, vector<double> z_in, vector<double> alpha_in);		//��ԭʼ�������룬��output�ɶ�
	vector<vector<double>> output();
    vector<vector<double>> outputByKmeans();																													//output ��input�ɶ�
	void set_dflist(vector<DFdata> dflist_in); // ����dflist����ʼ��
	void set_Pe(double pe_in);	//�������Ŷ�
	double get_Pe();			//��ȡ��ǰ���Ŷ�
	int get_n();			//��ȡ ��ǰ�������ݵ�����
	vector<double> get_locxyz(); // ��ȡ����xyz		  
	vector<double> get_ellipse(); //��ȡ�����Բ	
	vector<DFdata> get_dflist();	//��ȡdflist
	bool isValid();				//�Ƿ���Ч
	double AverageDistance;		//ƽ������
};
#endif