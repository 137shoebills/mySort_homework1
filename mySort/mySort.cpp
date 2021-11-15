#include <iostream>
#include <ctime>
#include <algorithm>
#include <fstream>

const int maxn = 400010;
double recur_num = 0;//递归层数全局变量
int ave_num = 0;//用来partition的平均值
int t[maxn] = {}, mx;//用来计算dd值的树状数组
/*
* void add(int x, int v)
* 树状数组，插入函数
*/
void add(int x, int v) {
    while (x <= mx) {
        t[x] += v;
        x += (x & -x);
    }
}
/*
* int getit(int x)
* 树状数组，读取函数
*/
int getit(int x) {
    int num = 0;
    while (x) {
        num += t[x];
        x -= (x & -x);
    }return num;
}
/*
* bool cmp(int x)
* 构造序列的partition所用的cmp
*/
bool cmp(int x) {
    return x < ave_num;
}
/*
* bool judge_arr(int arr[], int siz)
* arr为输入的数组的头指针，siz为数组长度
* 当输入数组arr[]已经按照非递减序排列好时，不作变更，返回1
* 当arr[]已经按照非递增序排列好时，将arr[]变为逆序，返回1
* 若arr无序，函数返回0
*/
bool judge_arr(int arr[], int left,int right) {
    int flag1 = 1, flag2 = 1;
    for (int i = left; i <= right; ++i) {
        if (arr[i] < arr[i - 1]) {
            flag1 = 0;
        }
        if (arr[i] > arr[i - 1]) {
            flag2 = 0;
        }
    }
    if (flag1) {
        return 1;
    }
    else if (flag2) {
        for (int i = left, j = right; i < j; ++i, --j) {
            int c = arr[i];
            arr[i] = arr[j];
            arr[j] = c;
        }
        return 1;
    }return 0;
}
/*
* void my_merge1(int arr[],int left,int right,int temp[], int k)
* arr为排序的数组的头指针，left为当前处理区间的左序号，right为当前处理区间的右序号，temp为合并排序所需的额外空间,k为递归层数
* 递归合并排序的递归排序部分
*/
void my_merge1(int arr[],int left,int right,int temp[], int k) {
    recur_num = k > recur_num ? k : recur_num;
    if (left >= right)return;
    if (judge_arr(arr, left, right))return;
    int mid = (left + right) / 2;
    my_merge1(arr, left, mid, temp, k + 1);
    my_merge1(arr, mid + 1, right, temp, k + 1);
    int i = left, j = mid + 1, cnt = left;
    for (; cnt <= right;) {
        if (i <= mid && (j > right || arr[i] < arr[j])) {
            temp[cnt++] = arr[i];
            ++i;
        }
        else {
            temp[cnt++] = arr[j];
            ++j;
        }
    }
    for (i = left; i <= right; ++i) {
        arr[i] = temp[i];
    }
}
/*
* void my_merge2(int arr[],int length,int temp[])
* arr为排序的数组的头指针，length为数组长度，temp为合并排序所需的额外空间
* 非递归合并排序的排序部分
*/
void my_merge2(int arr[],int length,int temp[]) {
    for (int k = 2; k < length * 2; k *= 2) {
        int cnt = 0;
        for (int left = 0; left < length; left += k) {
            int right = left + k - 1;
            int mid = left + k / 2 - 1;
            if (mid >= length-1)break;
            if (right > length)right = length - 1;
            int i = left, j = mid + 1;
            for (; cnt <= right;) {
                if (i <= mid && (j > right || arr[i] < arr[j])) {
                    temp[cnt++] = arr[i];
                    ++i;
                }
                else {
                    temp[cnt++] = arr[j];
                    ++j;
                }
            }
        }
        for (int i = 0; i < length; ++i) {
            arr[i] = temp[i];
        }
    }
}
/*
* double merge_sort1(int arr[], int siz, int b[])
* arr为排序的数组的头指针，siz为数组长度,b为多次测试用的复原数组
* 递归合并排序
* 返回值为排序所用时间
*/
double merge_sort1(int arr[], int siz, int b[]) {
    clock_t start = clock();
    int* temp = new int[siz];
    double add_num = 0;
    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i < siz; ++i)arr[i] = b[i];
        recur_num = 0;
        my_merge1(arr, 0, siz - 1, temp, 1);
        add_num += recur_num;
    }recur_num = add_num/100.0;
    delete[]temp;
    return double(clock() - start) / 100.0;
}
/*
* double merge_sort2(int arr[], int siz, int b[])
* arr为排序的数组的头指针，siz为数组长度,b为多次测试用的复原数组
* 非递归合并排序
* 返回值为排序所用时间
*/
double merge_sort2(int arr[], int siz, int b[]) {
    clock_t start = clock();
    int* temp = new int[siz];
    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i < siz; ++i)arr[i] = b[i];
        recur_num = 0;
        my_merge2(arr, siz, temp);
    }
    delete[]temp;
    return double(clock() - start)/100.0;
}
/*
* void my_quick1(int arr[], int left, int right, int k)
* arr为排序的数组的头指针，left为当前处理区间的左序号，right为当前处理区间的右序号,k为递归层数
* 第一种快速排序的递归排序部分
*/
void my_quick1(int arr[], int left, int right, int k) {
    recur_num = k > recur_num ? k : recur_num;
    if (left >= right)return;
    if (judge_arr(arr, left, right))return;
    int l=left, h=right;
    int num = arr[left];
    for (; l < h;) {
        while (l < h && arr[h] >= num) {
            --h;
        }arr[l] = arr[h];
        while (l < h && arr[l] <= num) {
            ++l;
        }arr[h] = arr[l];
    }arr[l] = num;
    my_quick1(arr, left, l-1, k + 1);
    my_quick1(arr, l+1, right, k + 1);
}
/*
* void my_quick2(int arr[], int left, int right, int k)
* arr为排序的数组的头指针，left为当前处理区间的左序号，right为当前处理区间的右序号,k为递归层数
* 第二种快速排序的递归排序部分
*/
void my_quick2(int arr[], int left, int right, int k) {
    recur_num = k > recur_num ? k : recur_num;
    if (left >= right)return;
    if (judge_arr(arr, left, right))return;
    int l = left, h = right;
    int id = left + rand() % (right - left + 1);
    int num = arr[id]; 
    arr[id] = arr[left];
    arr[left] = num;
    for (; l < h;) {
        while (l < h && arr[h] >= num) {
            --h;
        }arr[l] = arr[h];
        while (l < h && arr[l] <= num) {
            ++l;
        }arr[h] = arr[l];
    }arr[l] = num;
    my_quick2(arr, left, l - 1, k + 1);
    my_quick2(arr, l + 1, right, k + 1);
}
/*
* double quick_sort1(int arr[], int siz , int b[])
* arr为排序的数组的头指针，siz为数组长度,b为多次测试用的复原数组
* 划分基准元素x固定选取当前区间最左端元素的快速排序1
* 返回值为排序所用时间
*/
double quick_sort1(int arr[], int siz , int b[]) {
    clock_t start = clock();
    double add_num = 0;
    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i < siz; ++i)arr[i] = b[i];
        recur_num = 0;
        my_quick1(arr, 0, siz - 1, 1);
        add_num += recur_num;
    }recur_num = add_num / 100.0;
    return double(clock() - start) / 100.0;
}
/*
* double quick_sort2(int arr[], int siz, int b[])
* arr为排序的数组的头指针，siz为数组长度,b为多次测试用的复原数组
* 划分基准元素x随机选取当前区间元素的快速排序2
* 返回值为排序所用时间
*/
double quick_sort2(int arr[], int siz, int b[]) {
    clock_t start=clock();
    srand((unsigned int)time(NULL));
    double add_num = 0;
    for (int k = 0; k < 100; ++k) {
        for (int i = 0; i < siz; ++i)arr[i] = b[i];
        recur_num = 0;
        my_quick2(arr, 0, siz - 1, 1);
        add_num += recur_num;
    }recur_num = add_num / 100.0;
    return double(clock() - start) / 100.0;
}
/*
* long long getDD(int arr[], int length) 
* arr为排序的数组的头指针，length为数组长度
* 通过树状数组在O(NlogN)的复杂度下计算序列DD值
* 返回值为序列DD值
*/
long long getDD(int arr[], int length) {
    mx = length;
    memset(t, 0, sizeof(t));
    long long num = 0;
    for (int i = length - 1; i >= 0; --i) {
        num += getit(arr[i]);
        add(arr[i]+1,1);
    }
    return num;
}
/*
* void getSequence(int arr[], int length, long long min_num, long long max_num, long long &dd) 
* arr为排序的数组的头指针，length为数组长度，min_num为dd值下限，max_num为上限，dd为初始复杂度
* 得到一个dd值在所求上限与下限之间的序列
*/
void getSequence(int arr[], int length, long long min_num, long long max_num, long long &dd) {
    int l = 0, ls = 0;
    int r = length - 1, rs = length - 1;
    int unit = length;
    while (dd < min_num || dd>max_num) {
        if(dd < min_num) {
            std::random_shuffle(arr + ls, arr + rs + 1);
            dd = getDD(arr, length);
            ls += unit, rs += unit;
            if (rs >= length)rs = length - 1;
            if (ls >= length) {
                if (unit == 1)break;
                ls = 0;
                rs = unit - 1;
                continue;
            }
        }
        if(dd > max_num) {
            ave_num = 0;
            for (int i = l; i <= r; ++i)ave_num += arr[i];
            ave_num /= (r - l + 1);
            std::partition(arr + l, arr + r + 1, cmp);
            dd = getDD(arr, length);
            l += unit, r += unit;
            if (r >= length)r = length - 1;
            if (l >= length) {
                unit /= 2;
                if (unit == 1)break;
                l = 0;
                r = unit - 1;
                continue;
            }
        }
    }
}
int main(){
    std::ofstream sequenceTXT;
    sequenceTXT.open("sequence.txt", std::ios::out);
    std::ofstream merge1TXT;
    merge1TXT.open("merge1.txt", std::ios::out);
    std::ofstream merge2TXT;
    merge2TXT.open("merge2.txt", std::ios::out);
    std::ofstream quick1TXT;
    quick1TXT.open("quick1.txt", std::ios::out);
    std::ofstream quick2TXT;
    quick2TXT.open("quick2.txt", std::ios::out);
    std::ofstream answerTXT;
    answerTXT.open("answer.txt", std::ios::out); 
    int len[6] = {2000,5000,10000,20000,30000,50000};
    for (int k = 0; k < 6; ++k) {
        long long length = len[k], maxdd, mindd;

        std::cout << len[k] << std::endl;
        answerTXT <<"正整数序列长度:"<<len[k]<<std::endl;
        int *a = new int[length+10];
        for (int i = 0; i < length; ++i)a[i] = i+1;
        std::random_shuffle(a, a + length);
        long long dd = getDD(a, length);
        answerTXT << "随机洗牌后序列的参考DD值为：" << dd << std::endl;
        maxdd = dd - 1; mindd = maxdd / 2;

        for (int w = 0; ; ++w) {
            if (double(dd) / double(length) < 50)break;
            if (w != 0) {
                getSequence(a, length, mindd, maxdd, dd);;
                maxdd = dd /2; mindd = maxdd / 2;
            }
            answerTXT << "当前序列DD值为：" << dd << std::endl;
            answerTXT << "当前序列ADD值为：" << double(dd)/double(length) << std::endl;
            for (int i = 0; i < length; ++i)sequenceTXT << a[i] << ' ';
            sequenceTXT << '\n';
            int* b = new int[length];

            answerTXT << "递归合并排序排序时间为" << merge_sort1(b, length, a) << std::endl;
            answerTXT << "递归层次为" << recur_num << std::endl;
            for (int i = 0; i < length; ++i)merge1TXT << b[i] << ' ';
            merge1TXT << '\n';

            answerTXT << "非递归合并排序排序时间为" << merge_sort2(b, length, a) << std::endl;
            for (int i = 0; i < length; ++i)merge2TXT << b[i] << ' ';
            merge2TXT << '\n';

            answerTXT << "快速排序2排序时间为" << quick_sort2(b, length, a) << std::endl;
            answerTXT << "递归层次为" << recur_num << std::endl;
            for (int i = 0; i < length; ++i)quick1TXT << b[i] << ' ';
            quick1TXT << '\n';

            answerTXT << "快速排序1排序时间为" << quick_sort1(b, length, a) << std::endl;
            answerTXT << "递归层次为" << recur_num << std::endl;
            for (int i = 0; i < length; ++i)quick1TXT << b[i] << ' ';
            quick1TXT << '\n';
            answerTXT << std::endl;
        }
        answerTXT  << std::endl;
    }
    return 0;
}
