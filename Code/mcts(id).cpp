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
#include<map>
using namespace std;
typedef unsigned int UI;
typedef long long LL;
typedef pair<int, int> PII;
typedef pair<UI, int>PUI;
typedef pair<int, pair<int, int> > PIII;
typedef pair<UI, UI> PUU;
#define yx second.first
#define yy second.second
#define x first
#define y second
const UI INF = 0xffffffffll;
struct TimeWindow
{
	UI ts, tc;
	int u;
	bool operator<(const TimeWindow& TW)const
	{
		if (ts != TW.ts)return ts < TW.ts;
		if (tc != TW.tc)return tc < TW.tc;
		return u < TW.u;
	}
};
struct MCTS_Link
{
	UI ts, tc;
	int ne;
};
vector<int> v;
int vern, kmax;
LL sum2;
UI tmax;
vector<vector<map<UI, UI> > > PHC_Index;
vector<int> core;
vector<vector<TimeWindow> >k_timewindow;
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
vector<int>transnode;
void build_MCTS_K(int k)
{
	transnode.resize(vern + 2);
	printf("start build MCTS-index:%d\n", k + 1);
	vector<UI>rectw(vern + 2, 0);
	vector<int>recnode;
	recnode.push_back(0);
	recnode.push_back(vern + 1);
	for (int i = 1; i <= vern; i++)
	{
		if (core[i] > k)
		{
			recnode.push_back(i);
		}
	}
	sort(recnode.begin(), recnode.end());
	for (int i = 0; i < recnode.size(); i++)
	{
		transnode[recnode[i]] = i;
	}
	MCTS_index[k].resize(recnode.size());
	set<PUI> twordertree;
	twordertree.insert({ 0,0 });
	twordertree.insert({ INF,vern + 1 });
	rectw[0] = 0;
	rectw[vern + 1] = INF;
	vector<TimeWindow>& k_tws = k_timewindow[k];
	MCTS_index[k][transnode[0]].push_back({ 0,INF,transnode[vern + 1] });
	sort(k_tws.begin(), k_tws.end());
	int n = k_tws.size();
	for (int i = 0; i < n; i++)
	{
		TimeWindow tw = k_tws[i];
		int u = tw.u;
		UI ts = tw.ts, tc = tw.tc;
		if (rectw[u] != 0)
		{
			set<PUI>::iterator it = twordertree.lower_bound({ rectw[u],u });
			set<PUI>::iterator preit = it;
			preit--;
			set<PUI>::iterator pastit = it;
			pastit++;
			int preu = (*preit).y;
			int pastu = (*pastit).y;
			UI pasttc = (*pastit).x;
			if (MCTS_index[k][transnode[preu]].size() == 0 || MCTS_index[k][transnode[preu]].back().ts != ts)
				MCTS_index[k][transnode[preu]].push_back({ ts,pasttc,transnode[pastu] });
			else
			{
				MCTS_index[k][transnode[preu]].back().tc = pasttc;
				MCTS_index[k][transnode[preu]].back().ne = transnode[pastu];
			}
			twordertree.erase(it);
			rectw[u] = tc;
			if (tc == INF)continue;
			PUI tc_u = make_pair(tc, u);
			it = twordertree.insert(tc_u).first;
			preit = it;
			preit--;
			pastit = it;
			pastit++;
			preu = (*preit).y;
			pastu = (*pastit).y;
			pasttc = (*pastit).x;
			if (MCTS_index[k][transnode[preu]].size() == 0 || MCTS_index[k][transnode[preu]].back().ts != ts)
				MCTS_index[k][transnode[preu]].push_back({ ts,tc,transnode[u] });
			else
			{
				MCTS_index[k][transnode[preu]].back().tc = tc;
				MCTS_index[k][transnode[preu]].back().ne = transnode[u];
			}
			if (MCTS_index[k][transnode[u]].size() == 0 || MCTS_index[k][transnode[u]].back().ts != ts)
				MCTS_index[k][transnode[u]].push_back({ ts,pasttc,transnode[pastu] });
			else
			{
				MCTS_index[k][transnode[u]].back().tc = pasttc;
				MCTS_index[k][transnode[u]].back().ne = transnode[pastu];
			}
		}
		else
		{
			rectw[u] = tc;
			PUI tc_u = make_pair(tc, u);
			set<PUI>::iterator it = twordertree.insert(tc_u).first;
			set<PUI>::iterator preit = it;
			preit--;
			set<PUI>::iterator pastit = it;
			pastit++;
			int preu = (*preit).y;
			int pastu = (*pastit).y;
			UI pasttc = (*pastit).x;
			if (MCTS_index[k][transnode[preu]].size() == 0 || MCTS_index[k][transnode[preu]].back().ts != ts)
				MCTS_index[k][transnode[preu]].push_back({ ts,tc,transnode[u] });
			else
			{
				MCTS_index[k][transnode[preu]].back().tc = tc;
				MCTS_index[k][transnode[preu]].back().ne = transnode[u];
			}
			if (MCTS_index[k][transnode[u]].size() == 0 || MCTS_index[k][transnode[u]].back().ts != ts)
				MCTS_index[k][transnode[u]].push_back({ ts,pasttc,transnode[pastu] });
			else
			{
				MCTS_index[k][transnode[u]].back().tc = pasttc;
				MCTS_index[k][transnode[u]].back().ne = transnode[pastu];
			}
		}
	}
}
void build_MCTS()
{
	printf("start build TimeWindowList\n");
	auto tstart = clock();
	k_timewindow.resize(kmax);
	MCTS_index.resize(kmax);
	for (int i = 1; i <= vern; i++)
	{
		for (int j = 0; j < core[i]; j++)
		{
			for (PUU tw : PHC_Index[i][j])
				k_timewindow[j].push_back({ tw.x,tw.y,i });
		}
	}
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
