#include"PLandEE.h"
using namespace std;
//数据转换
vector<double> GeographicToGeodeticRectangular(double Long, double Lati, double Height){
	double dN = 0;
	double EARTH_a = 6378137;                               // 地球赤道半径，单位：米
	double EARTH_e = 0.0818191908;                          // 地球第一偏心率                   // 地球第一偏心率
	double EARTH_b = EARTH_a * sqrt(1 - EARTH_e * EARTH_e); // 地球极半径，单位：米
	double EARTH_E_2 = pow(EARTH_e, 2);

	Lati = Lati * acos(-1) / 180;
	Long = Long * acos(-1) / 180;
	dN = EARTH_a / sqrt(1 - EARTH_E_2 * sin(Lati) * sin(Lati));
	double x = (dN + Height) * cos(Lati) * cos(Long);
	double y = (dN + Height) * cos(Lati) * sin(Long);
	double z = (dN * (1 - EARTH_E_2) + Height) * sin(Lati);
	vector<double> res;
	res.push_back(x); res.push_back(y); res.push_back(z);
	return res;
}
vector<double> GeodeticRectangularToStationRectangular(double xg, double yg, double zg, double Long, double Lati, double Height){
	double TransMatrix[3][3];
	double N;
	double EARTH_a = 6378137;                               // 地球赤道半径，单位：米
	double EARTH_e = 0.0818191908;                          // 地球第一偏心率                   // 地球第一偏心率
	double EARTH_b = EARTH_a * sqrt(1 - EARTH_e * EARTH_e); // 地球极半径，单位：米
	double EARTH_E_2 = pow(EARTH_e, 2);
	double EARTH_E1_2 = (pow(EARTH_a, 2) - pow(EARTH_b, 2)) / pow(EARTH_b, 2);
	double DEG2RAD = 0.017453292519943295769236907684886;
	double EPOWER = 0.0066943800042608065155834872766664;
	Long *= acos(-1) / 180;
	Lati *= acos(-1) / 180;
	N = EARTH_a / sqrt(1 - EARTH_E_2 * (pow(sin(Lati), 2)));//卯酉圈半径
	TransMatrix[0][0] = -sin(Long);
	TransMatrix[0][1] = cos(Long);
	TransMatrix[0][2] = 0;
	TransMatrix[1][0] = -cos(Long) * sin(Lati);
	TransMatrix[1][1] = -sin(Long) * sin(Lati);
	TransMatrix[1][2] = cos(Lati);
	TransMatrix[2][0] = cos(Long) * cos(Lati);
	TransMatrix[2][1] = sin(Long) * cos(Lati);
	TransMatrix[2][2] = sin(Lati);
	double xStation = TransMatrix[0][0] * (xg - (N + Height) * cos(Long) * cos(Lati)) +
		TransMatrix[0][1] * (yg - (N + Height) * sin(Long) * cos(Lati)) +
		TransMatrix[0][2] * (zg - (N * (1 - EARTH_E_2) + Height) * sin(Lati));
	double yStation = TransMatrix[1][0] * (xg - (N + Height) * cos(Long) * cos(Lati)) +
		TransMatrix[1][1] * (yg - (N + Height) * sin(Long) * cos(Lati)) +
		TransMatrix[1][2] * (zg - (N * (1 - EARTH_E_2) + Height) * sin(Lati));
	double zStation = TransMatrix[2][0] * (xg - (N + Height) * cos(Long) * cos(Lati)) + TransMatrix[2][1] * (yg - (N + Height) * sin(Long) * cos(Lati)) + TransMatrix[2][2] * (zg - (N * (1 - EARTH_E_2) + Height) * sin(Lati));
	vector<double> res;
	res.push_back(xStation); res.push_back(yStation); res.push_back(zStation);
	return res;
}
vector<double> StationRectangularToGeodeticRectangular(double xStationTarget, double yStationTarget, double zStationTarget, double longitudeStation, double latitudeStation, double heightStation){
	double N0;//,ePower,DEG2RAD;
	double x, y, z;

	//ePower  = 0.0066943800042608065155834872766664;
	//DEG2RAD = 0.017453292519943295769236907684886;
	double EARTH_a = 6378137;                               // 地球赤道半径，单位：米
	double EARTH_e = 0.0818191908;                          // 地球第一偏心率                   // 地球第一偏心率
	double EARTH_b = EARTH_a * sqrt(1 - EARTH_e * EARTH_e); // 地球极半径，单位：米
	double EARTH_E_2 = pow(EARTH_e, 2);
	double EARTH_E1_2 = (pow(EARTH_a, 2) - pow(EARTH_b, 2)) / pow(EARTH_b, 2);
	double DEG2RAD = 0.017453292519943295769236907684886;
	double EPOWER = 0.0066943800042608065155834872766664;
	longitudeStation = longitudeStation * DEG2RAD;
	latitudeStation = latitudeStation * DEG2RAD;
	N0 = EARTH_a / sqrt(1 - EPOWER * sin(latitudeStation) * sin(latitudeStation));//卯酉圈半径

	x = (N0 + heightStation) * cos(longitudeStation) * cos(latitudeStation);
	y = (N0 + heightStation) * sin(longitudeStation) * cos(latitudeStation);
	z = (N0 * (1 - EPOWER) + heightStation) * sin(latitudeStation);

	double X = -sin(longitudeStation) * xStationTarget - cos(longitudeStation) * sin(latitudeStation) * yStationTarget + cos(longitudeStation) * cos(latitudeStation) * zStationTarget + x;
	double Y = cos(longitudeStation) * xStationTarget - sin(longitudeStation) * sin(latitudeStation) * yStationTarget + sin(longitudeStation) * cos(latitudeStation) * zStationTarget + y;
	double Z = cos(latitudeStation) * yStationTarget + sin(latitudeStation) * zStationTarget + z;
	vector<double> res;
	res.push_back(X);
	res.push_back(Y);
	res.push_back(Z);
	return res;
}
vector<double> GeodeticRectangularToGeographic(double dX, double dY, double dZ){
	double p = 0;
	double seita = 0;;
	double EARTH_a = 6378137;                               // 地球赤道半径，单位：米
	double EARTH_e = 0.0818191908;                          // 地球第一偏心率                   // 地球第一偏心率
	double EARTH_b = EARTH_a * sqrt(1 - EARTH_e * EARTH_e); // 地球极半径，单位：米
	double EARTH_E_2 = pow(EARTH_e, 2);
	double EARTH_E1_2 = (pow(EARTH_a, 2) - pow(EARTH_b, 2)) / pow(EARTH_b, 2);
	double DEG2RAD = 0.017453292519943295769236907684886;
	double EPOWER = 0.0066943800042608065155834872766664;
	// 采用 WGS-84 标准
	p = sqrt(pow(dX, 2) + pow(dY, 2));
	seita = atan((dZ * EARTH_a) / (p * EARTH_b));

	double longitude = atan2(dY, dX);
	double latitude = atan((dZ + EARTH_E1_2 * EARTH_b * pow(sin(seita), 3)) / (p - EARTH_E_2 * EARTH_a * pow(cos(seita), 3)));
	double height = p / cos(latitude) - EARTH_a / sqrt(1 - EARTH_E_2 * pow(sin(latitude), 2));

	latitude = latitude * 180 / acos(-1);
	longitude = longitude * 180 / acos(-1);
	vector<double> res;
	res.push_back(longitude);
	res.push_back(latitude);
	res.push_back(height);
	return res;

}

double PicAngle2AlgAngle(double in)
{
	double res = 0;
	if (in <= 90){
		res = 90 - in;
	}
	else{
		in -= 90;
		res = 360 - in;
	}
	return res;
}
double AlgAngle2PicAngle(double in)
{
	double res = 0;
	if (in <= 90){
		res = 90 - in;
	}
	else{
		res = 450 - in;
	}
	return res;
}
/*
#####     #
#     #   ##
#     #  # #
#     #    #
#   # #    #
#    #     #
#### #  #####
*/
#define PI acos(-1)
//设置函数
void DFdata::set_x(double x_in){ x = x_in; }
void DFdata::set_y(double y_in){ y = y_in; }
void DFdata::set_z(double z_in){ z = z_in; }
void DFdata::set_alpha(double alpha_in, bool isPI){
	if (isPI){ alpha = alpha_in; }
	else{
		alpha = alpha_in / 360 * 2 * PI;  //转换成弧度
	}
	if (alpha > 2 * PI || alpha<0){ cout << "角度不合法，请重新设置角度,0--2PI(0-360)" << endl; }
}
void DFdata::set_beta(double beta_in, bool isPI){
	if (isPI){ beta = beta_in; }
	else{
		beta = beta_in / 360 * 2 * PI;
	}
	if (beta<-PI / 2 || beta>PI / 2){ cout << "角度不合法，请重新设置角度,-PI/2---PI/2" << endl; }
}
double DFdata::get_alpha(){ return alpha; }
double DFdata::get_beta(){ return beta; }
double DFdata::set_distance(double x_in, double y_in){
	distance = sqrt((x - x_in)*(x - x_in) + (y - y_in)*(y - y_in));
	return distance;
}
double DFdata::get_distance(){ return distance; }
bool DFdata::isValid(double x_in, double y_in){
	double n1[2];
	n1[0] = cos(alpha);
	n1[1] = sin(alpha);
	double n2[2];
	n2[0] = x_in - x;
	n2[1] = y_in - y;
	//向量点乘为正才合理
	bool flag = (n1[0] * n2[0] + n1[1] * n2[1])  >= 0;
	return flag;
}
vector<double> DFdata::get_point(){
	vector<double> p(3);
	p[0] = x;
	p[1] = y;
	p[2] = z;
	return p;
}
//DFdata ---------------------------- 算法中r参数的函数
vector<double> DFdata::set_r(double x_in, double y_in){
	r[0] = x - x_in;
	r[1] = y - y_in;
	vector<double> p(2);
	p[0] = r[0];
	p[1] = r[1];
	return p;
}
vector<double> DFdata::get_r(){
	vector<double> p(2);	// 为了对r值改变时不影响原始数据
	p[0] = r[0];
	p[1] = r[1];
	/*
	cout << "r" << endl;
	cout << p[0] << " " << p[1] << endl << endl;
	*/
	return p;
}
//DFdata ------------------------------sk相关函数
vector<double> DFdata::get_sk(){
	vector<double> p(2);
	p[0] = 1 * cos(alpha);
	p[1] = 1 * sin(alpha);
	///////test
/*	p[0] = 1 * sin(alpha);
	p[1] = 1 * cos(alpha);*/



	return p;
}
vector<double> DFdata::get_skt(){
	vector<double> p(2);
	double t = alpha + PI / 2;
	p[0] = cos(t);
	p[1] = sin(t);


	////////test
/*	p[0] = 1 * sin(t);
	p[1] = 1 * cos(t);
	*/
/*	cout << "getskt" << endl;
	cout << "alpha " << alpha*180/PI << " " << (alpha+PI/2)*180/PI<<endl;
	cout << p[0] << " " << p[1] << endl;
	cout << cos(alpha)*cos(t) + sin(alpha)*sin(t) << endl;
	cout << endl << endl;*/
	return p;
}
//DFdata -------------------------函数设置  构造函数
DFdata::DFdata(){}
DFdata::DFdata(double x_in, double y_in, double z_in, double alpha_in, double beta_in, bool isPI){
	x = x_in;
	y = y_in;
	z = z_in;
	if (isPI){ alpha = alpha_in; beta = beta_in; }
	else{
		alpha = alpha_in / 360 * 2 * PI;
		beta = beta_in / 360 * 2 * PI;
	}
	if (alpha > 2 * PI || alpha<0){ cout << "角度不合法，请重新设置角度alhpa,0--2PI(0-360)" << endl; }
	if (beta< -PI / 2 || beta>PI / 2){ cout << "角度不合法，请重新设置角度beta,-PI/2---PI/2" << endl; }
	get_sk();	//生成sk
	get_skt();	//生成skt向量
	set_r(0, 0);//以经纬度为0为默认值
}
//DFdata ---------------------------原始数据参数
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 基于测量数据的无源定位解法  Q1 算法一  PLocation  Passive Location
// PLocation 基础的参数设置和获取  初始化所用函数
PLocation::PLocation(){ Pe = 0.05; }
PLocation::PLocation(vector<DFdata> dflist_in){
	Valid = true;
	//选取第一个 输入的原始数据的id作为批次号；
	n = 0;
	Pe = 0.05;
	x = 0; y = 0; z = 0; rl = 0; rs = 0; fan = 0;
	int size = dflist_in.size();
	// 传入数据以及算法矩阵准备
	for (int i = 0; i < size; i++){
		dflist.push_back(dflist_in[i]);
		A.push_back(dflist_in[i].get_skt());
		b.push_back(help_b(dflist_in[i].get_skt(), dflist_in[i].get_r()));
		n++;
	}
	
	if (n == 0){ cout << "无数据 " << endl; return; } // 避免报错
	// 求得 定位点
	set_locxy();
	set_z();
	// 求得误差椭圆
	set_ellipse();

}
void PLocation::input(vector<double> x_in, vector<double> y_in, vector<double> z_in, vector<double> alpha_in){

	n = 0;
	if (x_in.size() == y_in.size()&&y_in.size() == z_in.size()&&z_in.size() == alpha_in.size()){
		dflist.clear();
		//坐标转换需要第一个测向数据的平台经纬度
		stationloc.clear();
		for (int i = 0; i < x_in.size(); i++){
			if (stationloc.size() == 0){
				stationloc.push_back(x_in.at(i));
				stationloc.push_back(y_in.at(i));
				stationloc.push_back(z_in.at(i));
			}
			//经纬度 转 地心坐标系
			vector<double> geoxyz = GeographicToGeodeticRectangular(x_in.at(i), y_in.at(i), 0.0);
			//地心 转 站心坐标系
			vector<double> staxyz = GeodeticRectangularToStationRectangular(geoxyz.at(0), geoxyz.at(1), geoxyz.at(2),
				stationloc.at(0), stationloc.at(1), stationloc.at(2));
			double alpha = PicAngle2AlgAngle(alpha_in.at(i));
			//测向数据输入格式 x y z  alpha--xy上偏角  beta--z上偏角  0/1 是否为弧度  1表示是
			DFdata tem(staxyz.at(0), staxyz.at(1), staxyz.at(2), alpha, 0.0, false);
			dflist.push_back(tem);
			n++;
		}
	}
	else{
		cout << "输入数据数目不对应 " << x_in.size() << " " << y_in.size() << " " << z_in.size() << " " << alpha_in.size() << endl;
		return;
	}
}
vector<vector<double>> PLocation::output(){
	vector<vector<double>> res;
	Valid = true;
	x = 0; y = 0; z = 0; rl = 0; rs = 0; fan = 0;
	if (n <= 1){
		cout << "数据量不够" << endl;
	}
	vector<DFdata> temdf = dflist;
	A.clear();
	b.clear();
	for (int i = 0; i < n; i++){
		DFdata tem = temdf.at(i);
		A.push_back(tem.get_skt());
		b.push_back(help_b(tem.get_skt(), tem.get_r()));
	}

	set_locxy();
	set_z();
	// 求得误差椭圆
	set_ellipse();
	//获取算法计算结果
	vector<double> locxyz = get_locxyz();
	vector<double> ell = get_ellipse();
	//站心转地心
	vector<double> sta2geo = StationRectangularToGeodeticRectangular(locxyz.at(0), locxyz.at(1), locxyz.at(2), stationloc.at(0), stationloc.at(1), stationloc.at(2));
	// 误差椭圆尝试转换 之后会改
	vector<double> ellxygeo = StationRectangularToGeodeticRectangular(ell.at(0), ell.at(1), locxyz.at(2), stationloc.at(0), stationloc.at(1), stationloc.at(2));
	vector<double> ellxy = GeodeticRectangularToGeographic(ellxygeo.at(0), ellxygeo.at(1), ellxygeo.at(2));
	//地心转经纬度
	locxyz = GeodeticRectangularToGeographic(sta2geo.at(0), sta2geo.at(1), sta2geo.at(2));
	ell[2] = ell[2] / acos(-1) * 180;
	ell[2] = PicAngle2AlgAngle(ell[2]) + 90;
	while (ell[2] > 180){
		ell[2] -= 180;
	}

	res.clear();
	res.push_back(locxyz);
	ell[0] = ellxy[0];
	ell[1] = ellxy[1];
	res.push_back(ell);
	result = res;
	return result;
}


vector<vector<double>> PLocation::outputByKmeans(){
    vector<vector<double>> res;
    Valid = true;
    x = 0; y = 0; z = 0; rl = 0; rs = 0; fan = 0;
    if (n <= 1){
        cout << "数据量不够" << endl;
    }
    vector<DFdata> temdf = dflist;
    A.clear();
    b.clear();

    int seq_len = 20; //每20个点一组
    int cnt = 0;

    vector<Point> kmeansPoint;
    for (int i = 0; i < n; i++) {
        DFdata tem = temdf.at(i);
        A.push_back(tem.get_skt());
        b.push_back(help_b(tem.get_skt(), tem.get_r()));
        ++cnt;
        if ((cnt == seq_len) || (n < seq_len)) {
            auto loc = set_locxy();
            Point point(loc);
            kmeansPoint.push_back(point);
            A.clear();
            b.clear();
            cnt = 0;
        }
    }
    auto points = KMeans::start(kmeansPoint, 1); //调用kmeans算法
    x = points[0].data_[0], y = points[0].data_[1];
    set_z();
    // 求得误差椭圆
    set_ellipse();
    //获取算法计算结果
    vector<double> locxyz = get_locxyz();
    vector<double> ell = get_ellipse();
    //站心转地心
    vector<double> sta2geo = StationRectangularToGeodeticRectangular(locxyz.at(0), locxyz.at(1), locxyz.at(2), stationloc.at(0), stationloc.at(1), stationloc.at(2));
    // 误差椭圆尝试转换 之后会改
    vector<double> ellxygeo = StationRectangularToGeodeticRectangular(ell.at(0), ell.at(1), locxyz.at(2), stationloc.at(0), stationloc.at(1), stationloc.at(2));
    vector<double> ellxy = GeodeticRectangularToGeographic(ellxygeo.at(0), ellxygeo.at(1), ellxygeo.at(2));
    //地心转经纬度
    locxyz = GeodeticRectangularToGeographic(sta2geo.at(0), sta2geo.at(1), sta2geo.at(2));
    ell[2] = ell[2] / acos(-1) * 180;
    ell[2] = PicAngle2AlgAngle(ell[2]) + 90;
    while (ell[2] > 180){
        ell[2] -= 180;
    }

    res.clear();
    res.push_back(locxyz);
    ell[0] = ellxy[0];
    ell[1] = ellxy[1];
    res.push_back(ell);
    result = res;
    return result;
}

void PLocation::set_dflist(vector<DFdata> dflist_in){
	//选取第一个 输入的原始数据的id作为批次号；
	dflist.clear();
	A.clear();
	b.clear();
	Valid = true;
	n = 0; x = 0; y = 0; z = 0; rl = 0; rs = 0; fan = 0;
	int size = dflist_in.size();
	for (int i = 0; i < size; i++){
		dflist.push_back(dflist_in[i]);
		A.push_back(dflist_in[i].get_skt());
		b.push_back(help_b(dflist_in[i].get_skt(), dflist_in[i].get_r()));
		n++;
	}
/*	cout << "a size" << size << endl;
	cout << "A" << endl;
	for (int i = 0; i < A.size(); i++){
		for (int j = 0; j < A.at(0).size(); j++){
			cout << " " << A.at(i).at(j);
		}
		cout << endl;
	}
	cout << "b" << endl;
	for (int i = 0; i < b.size(); i++){
		cout << " " << b.at(i);
	}
	cout << endl;;*/
	if (n == 0){ cout << "无数据 " << endl; }
	set_locxy();
	set_z();
	// 求得误差椭圆
	set_ellipse();
}
// 获取参数 ---------------- debug用
vector<vector<double>> PLocation::get_A(){
	vector<vector<double>> tem = A;
	return tem;
}
vector<double> PLocation::get_b(){
	vector<double> tem = b;
	return tem;
}
int PLocation::get_n(){ return n; }
void PLocation::set_Pe(double pe_in){ Pe = 1 - pe_in; }
double PLocation::get_Pe(){ return 1 - Pe; }
vector<DFdata> PLocation::get_dflist(){ return dflist; }

// 获取定位点
vector<double> PLocation::get_locxyz(){
	vector<double> point;
	point.push_back(x); point.push_back(y); point.push_back(z);
	return point;
}
//获取误差椭圆 相关函数
vector<double> PLocation::get_ellipse(){
	vector<double> tem;
	tem.push_back(rl); tem.push_back(rs); tem.push_back(fan);
	return tem;
}


//------------------------------------- 算法实现辅助函数  private  基础参数
double PLocation::help_b(vector<double> skt, vector<double> r){
	double tem = skt[0] * r[0] + skt[1] * r[1];
/*	cout << "skt" << endl;
	cout << skt[0] << " " << skt[1] << endl;
	cout << "r" << endl;
	cout << r[0] << " " << r[1] << endl;
	cout << "b" << tem << endl; */
	return tem;
}

//获取 定位点
vector<vector<double>> PLocation::get_ATA(){
	vector<vector<double>> ATA;
	int len = A.size();
	for (int i = 0; i < 2; i++)
	{
		vector<double> tem;
		for (int j = 0; j < 2; j++)
		{
			double sum = 0;
			for (int k = 0; k < len; k++)
			{
				sum += A[k][i] * A[k][j];
			}
			tem.push_back(sum);
		}
		ATA.push_back(tem);
	}
/*	cout << "AT A" << endl;
	for (int i = 0; i < ATA.size(); i++){
		for (int j = 0; j < ATA.at(0).size(); j++){
			cout << " " << ATA.at(i).at(j);
		}
		cout << endl;
	}*/
	return ATA;
}
bool PLocation::ATAsuit(vector<vector<double>> ATA){
	// 判断ATA为 一个2*2矩阵
	int len = ATA.size();
	if (len != 2){
		return 0;
	}
	for (int i = 0; i < 2; i++){
		int temlen = ATA[i].size();
		if (temlen != 2){ return 0; }
	}
	if (ATA[0][0] * ATA[1][1] == ATA[0][1] * ATA[1][0]){
		return 0;
	}
	return 1;
}
vector<vector<double>> PLocation::get_reATA(vector<vector<double>> ATA){
	vector<vector<double>> reATA;
	if (!ATAsuit(ATA)){
		cout << "错误，数据无法得出结论，ATA不可逆或是不是2*2矩阵" << endl;
		return reATA;
	}
	vector<double> tem1;
	// 方便计算，将 ATA的4个参数 声明为 a b c d,  reATA的4个参数声明为 x11 x12 x21 x22
	// 通过计算得到 x11 = d/(ad - bc)  x21 = -c/(ad - bc)  x12 = -b/(ad - bc) x22 = a/(ad - bc)
	double a = ATA[0][0]; double b = ATA[0][1]; double c = ATA[1][0]; double d = ATA[1][1];
	/*cout << "ata" << endl;
	cout << "a" << a << " b" << b << endl;
	cout << "c" << c << " d" << d << endl;*/
	double x11 = d / (a*d - b*c);	double x21 = -c / (a*d - b*c);
	double x12 = -b / (a*d - b*c);	double x22 = a / (a*d - b*c);
	tem1.push_back(x11);	tem1.push_back(x12);	reATA.push_back(tem1);
	vector<double> tem2;
	tem2.push_back(x21);	tem2.push_back(x22);	reATA.push_back(tem2);
	
	/*cout << "reATA" << endl;
	for (int i = 0; i < reATA.size(); i++){
		for (int j = 0; j < reATA.at(0).size(); j++){
			cout << " " << reATA.at(i).at(j);
		}
		cout << endl;
	}*/

	return reATA;
}
vector<vector<double>> PLocation::get_pinvA(){
	vector<vector<double>> ATA = get_ATA();
	if (ATAsuit(ATA)){
		vector<vector<double>> reATA = get_reATA(ATA);
		vector<vector<double>> temA = A;
		vector<vector<double>> res;
		int len = temA.size();
		for (int i = 0; i < 2; i++)
		{
			vector<double> tem;
			for (int j = 0; j < len; j++)
			{
				double num = reATA[i][0] * A[j][0] + reATA[i][1] * A[j][1];
				tem.push_back(num);
			}
			res.push_back(tem);
		}
		pinvA = res;
	
/*		cout << "pinvA"<<len<<endl;
		for (int i = 0; i < pinvA.size(); i++){
			for (int j = 0; j < pinvA.at(0).size(); j++){
				cout << " " << pinvA.at(i).at(j);
			}
			cout << endl;
		}	*/

		return res;
	}
	else{
		vector<vector<double>> err;
		return err;
	}

}
vector<double> PLocation::set_locxy(){
	pinvA = get_pinvA();
	vector<double> loc;
	if (pinvA.size() == 0){
		Valid = false;
		x = -1;
		y = -1;
		return loc;
	}
	for (int i = 0; i < 2; i++){
		double tem = 0;
		for (int j = 0; j < n; j++){
			tem += pinvA[i][j] * b[j];
//			cout << "tem " << tem << " pinvA " << pinvA[i][j] << " b " << b[j];
		}
//		cout << " x y" << endl;
//		cout << tem << endl;
		loc.push_back(tem);
	}
	x = loc[0];	y = loc[1];
	return loc;
}

//获取误差椭圆 相关 函数
void PLocation::set_ellipse(){
	if (!Valid){
		return;
	}
	vector<double> dis2; //distance 的 平方  获取d^2矩阵
	AverageDistance = 0;
	for (int i = 0; i < n; i++){
		//估计值到观察占距离
		double dis = dflist[i].set_distance(x, y);
		dis2.push_back(dis*dis);
		AverageDistance += dis;
	}
	AverageDistance /= n;
	// pinvA * d * pinvA^T
	// sigma平方矩阵  方差矩阵 这里我们减少运算 增加代码

	sig11 = 0;
	sig12 = 0;//sig12 = sig21;
	sig22 = 0;

	double angle_rr = 2.0;
	double scalee = tan(angle_rr / 180 * acos(-1))* tan(angle_rr / 180 * acos(-1));
	
	for (int i = 0; i < n; i++){
		sig11 += pinvA[0][i] * dis2[i] * pinvA[0][i] * scalee;
		sig12 += pinvA[0][i] * dis2[i] * pinvA[1][i] * scalee;
		sig22 += pinvA[1][i] * dis2[i] * pinvA[1][i] * scalee;
	}
	/*cout << "sig11 " << sig11 << endl
		<< "sig12 " << sig12 << endl
		<< "sig22 " << sig22 << endl;*/
	double lam1 = 1.0 / 2 * (sig11 + sig22 + sqrt((sig11 - sig22)*(sig11 - sig22) + 4 * sig12 * sig12));
	double lam2 = 1.0 / 2 * (sig11 + sig22 - sqrt((sig11 - sig22)*(sig11 - sig22) + 4 * sig12 * sig12));
	
	rl = sqrt(-2 * lam1*log(1 - Pe));
	rs = sqrt(-2 * lam2*log(1 - Pe)); 
/*	cout << "rl  " << rl << endl
		<< "rs " << rs << endl;*/
//	fan = 1.0 / 2.0 * atan(2 * sig12 / (sig11 - sig22));
	fan = func_fan(sig11, sig12, sig22);
}
double PLocation::set_z(){
	z = 0;
	if (!Valid){
		return z;
	}
	double zmax = -1.79769313486231570E+308;
	double zmin = 1.79769313486231570E+308;


	for (int i = 0; i < n; i++){
		double distance = dflist[i].set_distance(x, y);
		vector<double> point = dflist[i].get_point();
		// 测量点和定位点的水平差
		double beta = dflist[i].get_beta();
		double k = tan(beta)*distance + point[2];//这里beta是 示向角 在-PI/2 -- PI/2之间
		zmax = k > zmax ? k : zmax;
		zmin = k < zmin ? k : zmin;
		z += k;
	}
	z = z / n;
	return z;
}
double PLocation::func_fan(double sig11, double sig12, double sig22){
	//论文方法
	
/*	double thisfan = atan2(2 * sig12, sig11 - sig22);
	thisfan /= 2;
	//fan = -fan;
	if (sig22 < sig11){
		thisfan += acos(-1) / 2;
	}
	*/
	
	double thisfan = 0;
	double a;
	a = 1.0;
	double testa = 1;
	double b = -sig11 - sig22;
	double c = sig11*sig22 - sig12*sig12;

	double delta = sqrt(b*b - 4 * a*c);
	double sigee = (-b + delta) / 2 * a;
	thisfan = atan2(sigee-sig11, sig12);
	if (sig22 > sig11){
		thisfan += acos(-1) / 2;
	}	
/*	cout << "sig11 " << sig11 << endl
		<< "sig12 " << sig12 << endl
		<< "sig22 " << sig22 << endl
		<< "a " << a << endl
		<< "b " << b << endl
		<< "c " << c << endl
		<< "delta " << delta << endl
		<< "sigee " << sigee << endl
		<< "atan2(1, 1)"<<atan2(1,1)<<endl
		<< "fan "<<thisfan*180/acos(-1)<<endl;*/
	return thisfan;
}
//判断是否为射线得到的点
bool PLocation::isValid(){
	if (!Valid){
		return false;
	}
	bool flag = true;
	for (int i = 0; i < dflist.size(); i++){
		if (!dflist.at(i).isValid(x, y)){
			flag = false;
			break;
		}
	}
	return flag;
}