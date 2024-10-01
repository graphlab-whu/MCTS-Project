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
int ans1[10010 * 25], ans2[10010 * 25], qsi[25], qk[25];
vector<int> v;
int vern, kmax;
UI tmax;
vector<vector<map<UI, UI> > > PHC_Index;
vector<int> core;
vector<vector<vector<MCTS_Link> > > MCTS_index;

struct arc {
	int src, dst;
	LL t;
};
vector<arc>arcs(EMAX);
int arcn = 0;
vector<int> hs;
vector<int> hd;
vector<int> ht;
vector<int> spre, snxt;
vector<int>dpre, dnxt;
vector<int> tpre, tnxt;
vector<int> telarc;
int idx = 0;

int head = -1, tail = -1;
vector<int> tsarc;
vector<int> tsnxt, tspre;
int idt = 0;

template<typename T>
void hash_combine(size_t& seed, T const& v)
{
	seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash
{
	template<typename T1, typename T2>
	size_t operator()(std::pair<T1, T2>const& p) const
	{
		size_t seed = 0;
		hash_combine(seed, p.first);
		hash_combine(seed, p.second);
		return seed;
	}
};
unordered_map<int, int> Mv;
unordered_map<pair<int, int>, int, pair_hash> Mc;
set<pair<int, int>> Hv;

void _init_arr(vector<int>& vec, int  size)
{
	vec.resize(size);
}
void initmem()
{
	_init_arr(hs, VMAX);
	_init_arr(hd, VMAX);
	_init_arr(ht, TMAX);
	_init_arr(spre, EMAX), _init_arr(snxt, EMAX);
	_init_arr(dpre, EMAX), _init_arr(dnxt, EMAX);
	_init_arr(tpre, EMAX), _init_arr(tnxt, EMAX);
	_init_arr(telarc, EMAX);

	_init_arr(tsarc, EMAX);
	_init_arr(tsnxt, EMAX);
	_init_arr(tspre, EMAX);

}
vector<UI>ts;
void loadgraph(const char* name)
{
	ifstream fin(name, ios::in);
	if (fin.is_open() == false) { printf("open graph %s fail\n", name); exit(1); }
	vector<int> v;
	LL tmin = 1e12;

	string l;
	stringstream sstream;
	while (getline(fin, l))
	{
		sstream.clear();
		sstream.str(l);
		int ut, vt;
		LL tt;
		sstream >> ut >> vt >> tt;
		if (ut == vt)continue;
		v.push_back(ut);
		v.push_back(vt);
		tmin = min(tmin, tt);
		arcs[arcn++] = { ut,vt,tt };
	}
	fin.close();

	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	vern = v.size();

	auto get = [&](int k) {
		return (lower_bound(v.begin(), v.end(), k) - v.begin()) + 1;
		};
	for (int i = 0; i < arcn; ++i)
	{
		arcs[i].t -= tmin - 1;
		ts.push_back(arcs[i].t);
	}
	sort(ts.begin(), ts.end());
	ts.erase(unique(ts.begin(), ts.end()), ts.end());
	for (int i = 0; i < arcn; ++i)
	{
		arcs[i].src = get(arcs[i].src), arcs[i].dst = get(arcs[i].dst);
		if (arcs[i].src > arcs[i].dst) swap(arcs[i].src, arcs[i].dst);
		arcs[i].t = lower_bound(ts.begin(), ts.end(), arcs[i].t) - ts.begin() + 1;
	}
}
void addarc(int id, int src, int dst, int t)
{
	telarc[idx] = id;
	snxt[idx] = hs[src]; if (hs[src] >= 0) spre[hs[src]] = idx; hs[src] = idx;
	dnxt[idx] = hd[dst]; if (hd[dst] >= 0) dpre[hd[dst]] = idx; hd[dst] = idx;
	tnxt[idx] = ht[t]; if (ht[t] >= 0) tpre[ht[t]] = idx; ht[t] = idx;
	idx++;
}
void addt(UI t)
{
	if (head == -1) head = idt;
	tsarc[idt] = t;
	tspre[idt] = tail, tsnxt[idt] = -1;
	if (tail >= 0) tsnxt[tail] = idt; tail = idt;
	idt++;
}

void buildtel(int l, int r)
{
	fill(hs.begin(), hs.end(), -1);
	fill(hd.begin(), hd.end(), -1);
	fill(ht.begin(), ht.end(), -1);
	head = -1, tail = -1;
	idx = 0, idt = 0;
	int ts = 0, cnt = 0;
	vector<UI> tts;
	for (int i = 0; i < arcn; ++i)
		if (arcs[i].t >= l && arcs[i].t <= r)
		{
			addarc(i, arcs[i].src, arcs[i].dst, arcs[i].t);
			tts.push_back(arcs[i].t);
		}
	sort(tts.begin(), tts.end());
	tts.erase(unique(tts.begin(), tts.end()), tts.end());
	for (auto t : tts) addt(t);
}
bool cAdd(int src, int dst)
{
	if (Mc.count({ src, dst }) == 0)
	{
		Mc[{src, dst}] = 1;
		return true;
	}
	Mc[{src, dst}]++;
	return false;
}

void vAdd(int v)
{
	if (Mv.count(v) == 0) Mv[v] = 0;
	Hv.erase({ Mv[v], v });
	Mv[v]++;
	Hv.insert({ Mv[v], v });
}

void initMH(int l, int r)
{
	Mv.clear();
	Mc.clear();
	Hv.clear();
	for (int i = 0; i < arcn; ++i)
	{
		int src = arcs[i].src;
		int dst = arcs[i].dst;
		int t = arcs[i].t;
		if (t < l || t > r) continue;
		if (cAdd(src, dst))
		{
			vAdd(src);
			vAdd(dst);
		}
	}
}
void delt(int t)
{
	int i = (lower_bound(tsarc.begin(), tsarc.begin() + idt, t) - tsarc.begin());
	if (head == i) head = tsnxt[i];
	else if (tsnxt[i] == -1) tsnxt[tspre[i]] = -1, tail = tspre[i];
	else {
		tsnxt[tspre[i]] = tsnxt[i];
		tspre[tsnxt[i]] = tspre[i];
	}
}

void delarc_l(vector<int>& head, vector<int>& next, vector<int>& pre, int i, int idx)
{
	if (head[i] == idx) head[i] = next[idx];
	else if (next[idx] == -1) next[pre[idx]] = -1;
	else {
		next[pre[idx]] = next[idx];
		pre[next[idx]] = pre[idx];
	}
}

bool cUpd(int src, int dst)
{
#ifdef _DEBUF
	if (Mc.count({ src, dst }) == 0) { printf("cUpd:empty update\n"); exit(1); }
#endif
	bool ret = false;
	Mc[{src, dst}]--;
	if (Mc[{src, dst}] == 0) { Mc.erase({ src, dst }); ret = true; }
	return ret;
}

void vUpd(int v)
{
#ifdef _DEBUG 
	if (Mv.count(v) == 0) { printf("vUpd:empty update\n"); exit(1); }
#endif
	int d = Mv[v];
	Hv.erase({ d, v });
	Hv.insert({ d - 1, v });
	Mv[v]--;
}

void decomp(int k)
{
	while (Hv.size() && (Hv.begin()->first < k))
	{
		auto nv = *(Hv.begin());
		Hv.erase(Hv.begin());
		int n = nv.first, v = nv.second;
		Mv.erase(v);
		unordered_set<int> nbr;
		for (int i = hs[v]; ~i; i = snxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			cUpd(arcs[id].src, arcs[id].dst);
			int u = arcs[id].src == v ? arcs[id].dst : arcs[id].src;
			nbr.insert(u);
		}
		for (int i = hd[v]; ~i; i = dnxt[i])
		{
			int id = telarc[i];
			delarc_l(hs, snxt, spre, arcs[id].src, i);
			delarc_l(hd, dnxt, dpre, arcs[id].dst, i);
			delarc_l(ht, tnxt, tpre, arcs[id].t, i);
			if (ht[arcs[id].t] == -1) delt(arcs[id].t);
			cUpd(arcs[id].src, arcs[id].dst);
			int u = arcs[id].src == v ? arcs[id].dst : arcs[id].src;
			nbr.insert(u);
		}
		for (auto u : nbr) vUpd(u);
	}
}



int getL(UI ql)
{
	int l = 0, r = ts.size() - 1;
	while (l < r)
	{
		int mid = l + r >> 1;
		if (ts[mid] < ql)l = mid + 1;
		else r = mid;
	}
	return l + 1;
}
int getR(LL qr)
{
	int l = 0, r = ts.size() - 1;
	while (l < r)
	{
		int mid = l + r + 1 >> 1;
		if (ts[mid] > qr)r = mid - 1;
		else l = mid;
	}
	return r + 1;
}

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
vector<UI>cts;
void remove(UI ts, int u, int k)
{

	UI ct = cts[u];
	int vp = recpre[u], vn = recne[u];
	recne[vp] = vn;
	recpre[vn] = vp;
	if (MCTS_index[vp][k].back().ts == ts)
	{
		MCTS_Link& mcts_link = MCTS_index[vp][k].back();
		mcts_link.tc = cts[vn];
		mcts_link.ne = vn;
	}
	else
	{
		MCTS_index[vp][k].push_back({ ts,cts[vn],vn });
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
void update(UI ts, UI tc, int u, int k)
{
	remove(ts, u, k);
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
		if (MCTS_index[v][k].back().ts == ts)
		{
			MCTS_Link& mcts_link = MCTS_index[v][k].back();
			mcts_link.tc = cts[u];
			mcts_link.ne = u;
		}
		else
		{
			MCTS_index[v][k].push_back({ ts,cts[u],u });

		}
		if (MCTS_index[u][k].back().ts == ts)
		{
			MCTS_Link& mcts_link = MCTS_index[u][k].back();
			mcts_link.tc = cts[recne[u]];
			mcts_link.ne = recne[u];
		}
		else
		{
			MCTS_index[u][k].push_back({ ts,cts[recne[u]],recne[u] });
		}
	}
	else
	{
		auto it = recCT.lower_bound(tc);
		it--;
		PUI pi = *it;
		recCT[tc] = u;
		int pre = pi.y;
		recpre[recne[pre]] = u;
		recne[u] = recne[pre];
		recne[pre] = u;
		recpre[u] = pre;
		if (MCTS_index[pre][k].back().ts == ts)
		{
			MCTS_Link& mcts_link = MCTS_index[pre][k].back();
			mcts_link.tc = cts[u];
			mcts_link.ne = u;
		}
		else
		{
			MCTS_index[pre][k].push_back({ ts,cts[u],u });
		}
		if (MCTS_index[u][k].back().ts == ts)
		{
			MCTS_Link& mcts_link = MCTS_index[u][k].back();
			mcts_link.tc = cts[recne[u]];
			mcts_link.ne = recne[u];
		}
		else
		{
			MCTS_index[u][k].push_back({ ts,cts[recne[u]],recne[u] });
		}
	}
}
void build_MCTS_K(int k)
{
	printf("start build %d MCTS-Index\n", k + 1);
	rectc.clear();
	recCT.clear();

	recpre.resize(vern + 2, 0);
	recne.resize(vern + 2, 0);
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
	if (init_arr[0].y != 0 || init_arr[init_n - 1].y != vern + 1)
	{

		cout << "init_error 0\n";
		cout << init_arr[0].x << " " << init_arr[0].y << "\n" << init_arr[init_n - 1].x << " " << init_arr[init_n - 1].y << "\n";
		exit(0);
	}
	for (int i = 1; i < init_n - 1; i++)
	{
		if (init_arr[i].x != init_arr[i + 1].x)
		{
			recCT[init_arr[i].x] = init_arr[i].y;
		}
		recpre[init_arr[i].y] = init_arr[i - 1].y;
		recne[init_arr[i].y] = init_arr[i + 1].y;
		MCTS_index[init_arr[i].y][k].push_back({ 1,init_arr[i + 1].x,init_arr[i + 1].y });
	}
	MCTS_index[0][k].push_back({ 1,init_arr[1].x,init_arr[1].y });
	recne[0] = init_arr[1].y;
	recpre[vern + 1] = init_arr[init_n - 2].y;
	for (UI ts : rects)
	{
		for (PUI upd : rectc[ts])
		{
			update(ts, upd.x, upd.y, k);
		}
	}
}
void build_MCTS()
{
	printf("start build MCTS-Index\n");
	auto tstart = clock();
	MCTS_index.resize(vern + 1);
	cts.resize(vern + 2, 0);
	for (int i = 1; i <= vern; i++)
	{
		MCTS_index[i].resize(core[i]);
	}
	MCTS_index[0].resize(kmax);
	double MB = 1 << 20;
	for (int k = 1; k < kmax; k++)
	{
		build_MCTS_K(k);
	}
	printf("finish build MCTS-Index when %lf s\n", ((clock() - tstart) / (double)CLOCKS_PER_SEC));

	LL sum_mcts = 0;
	for (int i = 0; i <= vern; i++)
	{
		for (int j = 0; j < MCTS_index[i].size(); j++)
		{
			sum_mcts += MCTS_index[i][j].size() * 12ll;
		}
	}

	printf("MCTS-Index size is %lf MB\n", (double)sum_mcts / MB);
}
vector<int> query_By_PHC(UI ts, UI te, int k)
{
	vector<int>ret;
	for (int u = 1; u <= vern; u++)
	{
		if (core[u] >= k)
		{
			map<UI, UI>::iterator it = PHC_Index[u][k - 1].upper_bound(ts);
			it--;
			if ((*it).y <= te)
			{
				ret.push_back(u);
			}
		}
	}
	return ret;
}
vector<int> query_By_MCTS(UI ts, UI te, int k)
{
	vector<int>ret;
	int u = 0;
	while (true)
	{
		vector<MCTS_Link>& mcts_u_k = MCTS_index[u][k - 1];
		int l = 0, r = mcts_u_k.size() - 1;
		while (l < r)
		{
			int mid = l + r + 1 >> 1;
			if (mcts_u_k[mid].ts > ts)r = mid - 1;
			else l = mid;
		}
		MCTS_Link& mcts_link = mcts_u_k[l];
		if (mcts_link.tc <= te)
		{
			u = mcts_link.ne;
			ret.push_back(u);
		}
		else
		{
			break;
		}
	}
	return ret;
}
vector<int> query_By_TCD(int k)
{
	decomp(k);
	vector<int>ret;
	for (pair<int, int> tmp : Hv)
	{
		ret.push_back(tmp.second);
	}
	return ret;
}

UI getRand()
{
	return ((rand() << 16) | rand());
}
int base = 0;
void generateTest(double window_rate, double k_rate, int group_size)
{
	srand((unsigned)time(NULL));
	UI window_size = tmax * window_rate;
	int k_val = (int)(kmax * k_rate - 0.0001) + 1;
	UI range = tmax - window_size;
	bool flag = (k_rate == 0.1);
	for (int i = 0; i < group_size; i++)
	{
		if (flag)
		{
			UI ts = rand() % range + 1;
			UI te = ts + window_size;
			queries[i + base] = { ts,te };
		}
		else
		{
			queries[i + base] = queries[i + base - group_size];
		}
		
//		auto mcts_end = clock();
//		auto phc_start = clock();
//		vector<int>phc_result = query_By_PHC(queries[i + base].ts, queries[i + base].te, k_val);
//		auto phc_end = clock();
//		phc_cost += phc_end - phc_start;
//
//		auto mcts_start = clock();
//		vector<int>mcts_result = query_By_MCTS(queries[i + base].ts, queries[i + base].te, k_val);
//		auto mcts_end = clock();
//		mcts_cost += mcts_end - mcts_start;
//		
//		int l = getL(queries[i + base].ts), r = getR(queries[i + base].te);
//		buildtel(l, r);
//		initMH(l, r);
//		auto tcd_start = clock();
//		vector<int>tcd_result = query_By_TCD(k_val);
//		auto tcd_end = clock();
//		tcd_cost += tcd_end - tcd_start;
//		if (phc_result.size() != mcts_result.size() || mcts_result.size() != tcd_result.size())
//		{
//			printf("query error ts = %lld, te = %lld ,k = %d\n", (LL)queries[i + base].ts, (LL)queries[i + base].te, k_val);
//			exit(0);
//		}
	}
	auto phc_start = clock();
	for (int i = 0; i < group_size; i++)
	{
		vector<int>phc_result = query_By_PHC(queries[i + base].ts, queries[i + base].te, k_val);
		ans1[i + base] = phc_result.size();
	}
	auto phc_end = clock();
	double phc_cost = phc_end - phc_start;
	auto mcts_start = clock();
	for (int i = 0; i < group_size; i++)
	{
		vector<int>mcts_result = query_By_MCTS(queries[i + base].ts, queries[i + base].te, k_val);
		ans2[i + base] = mcts_result.size();
	}
	auto mcts_end = clock();
	double mcts_cost = mcts_end - mcts_start;
	
	
	double tcd_cost = 0;
	for(int i = 0;i < group_size;i++)
	{
		if(ans1[i + base] != ans2[i + base])
		{
			cout << "error\n";
			exit(0);
		}
	}
	base += group_size;
	phc_cost *= 1e6 / CLOCKS_PER_SEC / group_size;
	mcts_cost *= 1e6 / CLOCKS_PER_SEC / group_size;
	tcd_cost *= 1e6 / CLOCKS_PER_SEC / group_size;
	printf("window_rate is %.2lf, k_rate is %.2lf\ngroup_size is %d\nphc_time_cost is %.6lf us\n,mcts_time_cost is %.6lf us\n,tcd_time_cost is %.6lf us \n\n\n"
	,window_rate,k_rate, group_size, phc_cost, mcts_cost,tcd_cost);
}
int getGroupSize(const char* groupsize)
{
	int len = strlen(groupsize);
	int ret = 0;
	for (int i = 0; i < len; i++)
		ret = ret * 10 + groupsize[i] - '0';
	return ret;
}
int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("./mcts [graph] [phc] [groupsize]\n");
		return 1;
	}
	initmem();
	loadgraph(argv[1]);
	loadPHC(argv[2]);
	int groupsize = getGroupSize(argv[3]);
	build_MCTS();
	
	for (double window_rate = 0.1; window_rate < 1; window_rate += 0.2)
	{
		for (double k_rate = 0.1; k_rate < 1; k_rate += 0.2)
		{
			generateTest(window_rate, k_rate, groupsize);
		}
	}
	return 0;

}
