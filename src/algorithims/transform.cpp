//
// Created by ChenYL on 2022/8/3.
//
#include "transform.h"

using namespace std;

// ����ת���ɽǶ�
double hu2angle(double hu)
{
    return hu*180/PI;
}

// �Ƕ�ת���ɻ���
double angle2hu(double angle)
{
    return angle/180*PI;
}

// ��׼���Ƕ�
double to360(double angle)
{
    while (angle > 360) angle -= 360;
    while (angle < 0) angle += 360;
    return angle;
}

// �󽻵�  �Ƕȵ�λ����
vector<double> cross_point(double x1, double y1, double alpha1,
                           double x2, double y2, double alpha2){
    vector<double> res;
    // �� alpha������ 360���ڵ�
     alpha1 = to360(alpha1);
     alpha2 = to360(alpha2);
    alpha1 = angle2hu(alpha1);
    alpha2 = angle2hu(alpha2);
    //�������ƽ�� ���ش���
    if (abs(alpha1 - alpha2) < 10e-5){
        return res;
    }
    double l = ((y2 - y1) * cos(alpha2) - (x2 - x1) * sin(alpha2)) /
               (sin(alpha1 - alpha2));
    double x12 = x1 + l * cos(alpha1);
    double y12 = y1 + l * sin(alpha1);
    res.push_back(x12);
    res.push_back(y12);
    return res;
}
//������Ƕȣ��ο�����https://zhidao.baidu.com/question/1954601155876548588.html
double get_angle(double x1, double y1, double x2, double y2)
{
    double angle_res;
    double xx, yy;
    xx = x2 - x1;
    yy = y2 - y1;
    if (abs(xx) <= 10e-8){
        return acos(-1) / 2.0;
    }

    angle_res = atan(yy / xx);

    if (xx < 0.0){
        angle_res += 180.0;
    }
    else if (yy < 0.0) {
        angle_res += 360.0;
    }
    return angle_res;
}

bool is_line_point(double platx, double platy, double alpha, double tarx, double tary){
    double x = tarx - platx;
    double y = tary - platy;

    //��׼���Ƕ�
    alpha = to360(alpha);
    alpha = angle2hu(alpha);


    double anglex = cos(alpha);
    double angley = sin(alpha);

    if (x*anglex + y*angley > 0){
        return true;
    }
    else{
        return false;
    }
}

// ����ʱ��
struct tm create_time(int year, int month, int date, int hour, int min, int sec){
    struct tm t;
    if (year < 0){
        cout << "create_time err  " << " year " << year<<endl;
    }
    else if (month < 1 || month>12){
        cout << "create_time err  " << " month " << month << endl;
    }
    else if (date < 1 || date>31){
        cout << "create_time err  " << " date " << date << endl;
    }
    else if (hour < 0 || hour>23){
        cout << "create_time err  " << " hour " << hour << endl;
    }
    else if (min < 0 || min >59){
        cout << "create_time err  " << " min " << min << endl;
    }
    else if (sec < 0 || sec>59){
        cout << "create_time err  " << " sec " << sec << endl;
    }
    // ��Ҫʱ��ṹ������������ ��ݴ�1900��ʼ  �·�Ϊ0-11
    t.tm_year = year-1900; t.tm_mon = month-1; t.tm_mday = date;
    t.tm_hour = hour; t.tm_min = min; t.tm_sec = sec;
    return t;
}

string to_str_time(struct tm t){
    string s = to_string(t.tm_year+1900) + "-" +
               to_string(t.tm_mon+1) + "-" +
               to_string(t.tm_mday) + " " +
               to_string(t.tm_hour) + ":" +
               to_string(t.tm_min) + ":" +
               to_string(t.tm_sec);
    return s;
}

long sub_time(struct tm t_early, struct tm t_late){
    long res = LONG_MIN;
    res = difftime(mktime(&t_late), mktime(&t_early));
    return res;
}


bool same_point(double x1, double y1, double x2, double y2){
    if (abs(x1 - x2) < 10e-3 && abs(y1 - y2) < 10e-3){
        return true;
    }
    return false;
}


// �㵽ֱ�ߵľ���  �Ƕ�Ϊ �� ��λ
double getPointLineDistance(double x1, double y1, double alpha1, double x2, double y2){
    double alpha = angle2hu(to360(alpha1));
    //��
    double A = tan(alpha);
    double B = -1;
    double C = y1 - x1*tan(alpha);
    double shang = abs(A*x2 + B*y2 + C);
    double xia = sqrt(A*A + B*B);
    double distance = shang / xia;
    return distance;

}


// Map �ӿ� ��ʱ���Ӧע�͵�
double getTwoPointDistance(double x1, double y1, double x2, double y2){
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}