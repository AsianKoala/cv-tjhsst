#include <bits/stdc++.h>

using namespace std;

void printV(vector<int> v) {
    for(int i=0; i<v.size(); i++) {
        cout<<v[i]<<" ";
    }
    cout<<endl;
}

int main() {
    vector<int> v = {1,2,3,4};

    srand(time(NULL));
    for(int i=v.size()-1; i>0; i--) {
        int j = rand() % (i + 1);
        int temp = v[i];
        v[i] = v[j];
        v[j] = temp;
    }

    printV(v);
}