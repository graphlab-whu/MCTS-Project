#pragma GCC optimize(1)
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<set>
#include<algorithm>
#include<ctime>
#include<stdlib.h>
#include<unordered_map>
#include<unordered_set>
#include<map>
#include<cstring>
const int VMAX = 4000000;
const int EMAX = 83000000;
const int TMAX = 56650000;
const int GroupSize = 10000;
const int TCDGroupSize = 10;
using namespace std;
typedef unsigned int UI;
typedef long long LL;
typedef unsigned long long ULL;
typedef pair<int, int> PII;
typedef pair<UI, int>PUI;
typedef pair<UI, UI> PUU;
typedef pair<int, pair<int, int> > PIII;
#define yx second.first
#define yy second.second
#define x first
#define y second
const UI INF = 0xffffffffll;
struct MCTS_Link
{
	UI ts, tc;
	int ne;
};
struct Query
{
	UI ts, te;
};
Query queries[10010 * 25];
int ans1[10010 * 25], ans2[10010 * 25], ans3[10010 * 25];
vector<int> v;
int vern, kmax;
UI tmax;
vector<vector<map<UI, UI> > > PHC_Index;
vector<int> core;
vector<vector<vector<MCTS_Link> > > MCTS_index;

void loadPHC(const char* phc_file)
{
	kmax = 0;
	ifstream fin(phc_file, ios::in);
	if (!fin.is_open())
	{
		printf("fail to load file %s\n", phc_file);
		exit(1);
	}
	printf("start load PHC-index %s\n", phc_file);
	string s;
	getline(fin, s);
	stringstream sstream;
	sstream.str(s);
	sstream >> vern >> tmax;
	core.resize(vern + 1);
	PHC_Index.resize(vern + 1);
	v.resize(vern + 1);
	getline(fin, s);
	sstream.clear();
	sstream.str(s);
	for (int i = 1; i <= vern; i++)
	{
		sstream >> v[i];
	}
	for (int i = 1; i <= vern; i++)
	{
		getline(fin, s);
		sstream.clear();
		sstream.str(s);
		sstream >> core[i];
		kmax = max(kmax, core[i]);
		PHC_Index[i].resize(core[i]);
		for (int j = 0; j < core[i]; j++)
		{
			getline(fin, s);
			sstream.clear();
			sstream.str(s);
			int t;
			sstream >> t;
			for (int k = 0; k < t; k++)
			{
				UI x, y;
				getline(fin, s);
				sstream.clear();
				sstream.str(s);
				sstream >> x >> y;
				if (y == INF)y--;
				PHC_Index[i][j][x] = max(PHC_Index[i][j][x], y);
			}
		}
	}
	printf("finish load PHC_index\n");
	fin.close();
}
unordered_map<UI, vector<PUI>> rectc;
map<UI, int> recCT;
vector<int>recpre;
vector<int>recne;
vector<int>transnode;
vector<UI>cts;
void remove(UI ts, int u, int k, vector<vector<MCTS_Link>>& mcts_k)
{
	int MCTS_index = 0;
	UI ct = cts[u];
	int vp = recpre[u], vn = recne[u];
	recne[vp] = vn;
	recpre[vn] = vp;
	if (mcts_k[transnode[vp]].back().ts == ts)
	{
		MCTS_Link& mcts_link = mcts_k[transnode[vp]].back();
		mcts_link.tc = cts[vn];
		mcts_link.ne = transnode[vn];
	}
	else
	{
		mcts_k[transnode[vp]].push_back({ ts,cts[vn],transnode[vn] });
	}
	if (cts[vp] != ct && cts[vn] != ct)
	{
		recCT.erase(ct);
	}
	else if (cts[vn] != ct)
	{
		recCT[ct] = vp;
	}
}
void update(UI ts, UI tc, int u, int k, vector<vector<MCTS_Link>>& mcts_k)
{
	int MCTS_index = 0;
	remove(ts, u, k, mcts_k);
	cts[u] = tc;
	if (tc == INF - 1)return;
	if (recCT.count(tc))
	{
		int v = recCT[tc];
		recCT[tc] = u;
		recpre[u] = v;
		recne[u] = recne[v];
		recpre[recne[u]] = u;
		recne[v] = u;
		if (mcts_k[transnode[v]].back().ts == ts)
		{
			MCTS_Link& mcts_link = mcts_k[transnode[v]].back();
			mcts_link.tc = cts[u];
			mcts_link.ne = transnode[u];
		}
		else
		{
			mcts_k[transnode[v]].push_back({ ts,cts[u],transnode[u] });

		}
		if (mcts_k[transnode[u]].back().ts == ts)
		{
			MCTS_Link& mcts_link = mcts_k[transnode[u]].back();
			mcts_link.tc = cts[recne[u]];
			mcts_link.ne = transnode[recne[u]];
		}
		else
		{
			mcts_k[transnode[u]].push_back({ ts,cts[recne[u]],transnode[recne[u]] });
		}
	}
	else
	{
		auto it = recCT.upper_bound(tc);
		it--;
		PUI pi = *it;
		recCT[tc] = u;
		int pre = pi.y;
		recpre[recne[pre]] = u;
		recne[u] = recne[pre];
		recne[pre] = u;
		recpre[u] = pre;
		if (mcts_k[transnode[pre]].back().ts == ts)
		{
			MCTS_Link& mcts_link = mcts_k[transnode[pre]].back();
			mcts_link.tc = cts[u];
			mcts_link.ne = transnode[u];
		}
		else
		{
			mcts_k[transnode[pre]].push_back({ ts,cts[u],transnode[u] });
		}
		if (mcts_k[transnode[u]].back().ts == ts)
		{
			MCTS_Link& mcts_link = mcts_k[transnode[u]].back();
			mcts_link.tc = cts[recne[u]];
			mcts_link.ne = transnode[recne[u]];
		}
		else
		{
			mcts_k[transnode[u]].push_back({ ts,cts[recne[u]],transnode[recne[u]] });
		}
	}
}

void build_MCTS_K(int k)
{
	vector<vector<MCTS_Link>>& mcts_k = MCTS_index[k];
	int MCTS_index;
	printf("start build %d MCTS-Index\n", k + 1);
	rectc.clear();
	recCT.clear();

	recpre.resize(vern + 2);
	recne.resize(vern + 2);
	transnode.resize(vern + 2);
	vector<PUI>init_arr;
	init_arr.push_back({ 0,0 });
	init_arr.push_back({ INF,vern + 1 });
	set<UI> rects;

	for (int i = 1; i <= vern; i++)
	{
		if (core[i] < k + 1)
		{
			cts[i] = INF - 1;
			continue;
		}
		map<UI, UI>& phc = PHC_Index[i][k];
		map<UI, UI>::iterator it = phc.begin();
		init_arr.push_back({ (*it).y ,i });
		cts[i] = (*it).y;
		it++;
		while (it != phc.end())
		{
			rectc[(*it).x].push_back({ (*it).y,i });
			rects.insert((*it).x);
			it++;
		}
	}
	cts[0] = 0;
	cts[vern + 1] = INF;
	recCT[0] = 0;
	recCT[INF] = vern + 1;
	sort(init_arr.begin(), init_arr.end());
	int init_n = init_arr.size();
	mcts_k.resize(init_n);
	for (int i = 0; i < init_n; i++)
	{
		transnode[init_arr[i].y] = i;
	}
	for (int i = 1; i < init_n - 1; i++)
	{
		if (init_arr[i].x != init_arr[i + 1].x)
		{
			recCT[init_arr[i].x] = init_arr[i].y;
		}
		recpre[init_arr[i].y] = init_arr[i - 1].y;
		recne[init_arr[i].y] = init_arr[i + 1].y;
		mcts_k[transnode[init_arr[i].y]].push_back({ 1,init_arr[i + 1].x,transnode[init_arr[i + 1].y] });
	}
	mcts_k[transnode[init_arr[0].y]].push_back({ 1,init_arr[1].x,transnode[init_arr[1].y] });
	recne[0] = init_arr[1].y;
	recpre[vern + 1] = init_arr[init_n - 2].y;
	for (UI ts : rects)
	{
		for (PUI upd : rectc[ts])
		{
			update(ts, upd.x, upd.y, k, mcts_k);
		}
	}
}
void build_MCTS()
{
	printf("start build MCTS-Index\n");
	auto tstart = clock();
	MCTS_index.resize(kmax);
	cts.resize(vern + 2, 0);
	double MB = 1 << 20;
	for (int k = 1; k < kmax; k++)
	{
		build_MCTS_K(k);
	}
	printf("finish build MCTS-Index when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));

	LL sum_mcts = 0;
	for (int i = 1; i < kmax; i++)
	{
		for (int j = 0; j < MCTS_index[i].size(); j++)
		{
			sum_mcts += MCTS_index[i][j].size() * 12ll;
		}
	}

	printf("MCTS-Index size is %lf MB\n", (double)sum_mcts / MB);
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		printf("./mcts [phc]\n");
		return 1;
	}
	loadPHC(argv[1]);
	build_MCTS();
	return 0;
}
