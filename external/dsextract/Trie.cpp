#include "Trie.h"

Trie::Trie() {
	root = NULL;
	ncount = 0;
}
Trie::~Trie() {
	if(root != NULL)
		makeEmpty(root);
}

void Trie::insert(unsigned int v, vector<unsigned int> & outlink){
	if(root != NULL){
		if (outlink.size() > 0 && outlink[0] != root->label){
			// can't add to a trie a vertex outlink that do not share first element
			return;
		}
	}
	insert(v, outlink, root);
}

void Trie::printV( vector<unsigned int> v){
	for(unsigned int i=0; i<v.size(); i++)
		cout<<v[i]<<" ";
	cout<<endl;
}

void Trie::insert(unsigned int v, vector<unsigned int> &outlink, TrieNode * t){
	TrieNode *aux;
	TrieNode *ptr = NULL;
	int i=0;
	TrieNode *prevParent = NULL;
	for(vector<unsigned int>::iterator it=outlink.begin(); it != outlink.end(); it++){
		unsigned int posLabel = *it;
		aux = find(posLabel, t);
		// posLabel is the label of the  node and v is the vertex to add to set of node
		if(i == 0 && aux == NULL){
			insert(v, posLabel, t);
			ptr = t;
			t->depth = 1;
			prevParent = t;
			//cout<<"root "<<posLabel<<" v "<<v<<" t posLabel "<<t->label<<"\n";
			t=t->firstChild;
		} else if(i > 0 && aux == NULL && t == NULL){
			insert(v, posLabel, t);
			//cout<<"adding firstchild of "<<prevParent->label<<" firstChild is t->label "<<t->label<<" posLabel "<<posLabel<<" v "<<v<<"\n";
			if(prevParent == NULL){
				cout<<"prevParent NULL error \n";
				exit(1);
			}
			prevParent->nChildren++;
			t->parent = prevParent;
			prevParent->firstChild = t;
			t->depth = prevParent->depth + 1; 
			prevParent = t;
			t = t->firstChild;
		} else if(i > 0 && aux == NULL && t != NULL){
			//cout<<" adding in nextSibling of "<<t->label<<" v "<<v<<" poslabel "<<posLabel<<"\n";
			TrieNode *tpr = aux;	
			insert(v, posLabel, tpr);
			(t->parent)->nChildren++;
			tpr->depth = (t->parent)->depth + 1;
			tpr->parent = t->parent;
			//t->nextSibling = tpr;
			while(t->nextSibling != NULL)
				t = t->nextSibling;
			//cout<<" adding in nextSibling of "<<t->label<<" v "<<v<<" poslabel "<<posLabel<<"\n";
			//cout<<" t label "<<t->label<<endl;
			t->nextSibling = tpr;
			prevParent = tpr;
			t = tpr->firstChild;
		} else if( aux != NULL){
			// this only add a new field vertex to the node trie, it does not create a trienode
			//cout<<" node found label "<<aux->label<<" adding poslabel "<<posLabel<<" v "<<v<<"\n";
			insert(v, posLabel, aux);
			prevParent = aux;
			t = aux->firstChild;
		} 
		i++;
		
	}
	if (root == NULL)
		root = ptr;

	//cout<<" printTrie "<<endl;
	//printTrie();
}	
void Trie::insert(unsigned int v, unsigned int posLabel, TrieNode * & t){
	if(t == NULL){
		//cout<<"inserting posLabel "<<posLabel<<" v "<<v<<"\n";
		//cout<<" new "<<endl;
		ncount++;
		t = new TrieNode(v, posLabel,  NULL, NULL);
	} else {
		addToSet(v, posLabel, t);
	}
}

void Trie::addToSet(unsigned int v, unsigned int posLabel, TrieNode * & t){
	t->vertexSet.push_back(v);
}


TrieNode *Trie::find(unsigned int p, TrieNode * t){
	if( t == NULL){
		//cout<<" find t NULL"<<endl;
		return NULL;
	}
	if( p == t->label){
		//cout<<" find  p igual t label "<<p<<endl;
		return t;
	}
	else {
		return find(p,t->nextSibling);
	}
	
}

		
unsigned int Trie::findKey(unsigned int x){
	TrieNode* t = find (x, root);
	if(t == NULL)
		return 0;
	else
		return t->label;
}


void Trie::printTrie(){
	if(isEmpty())
		cout<<" Arbol vacio"<<endl;
	else
		printTrieRec(root);
}

void Trie::generateVirtualNodeList(vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN){
	vector<unsigned int> vnOut;
	vnOut.push_back(root->label);
	generateListIt(root,potentialVN,vnOut);
	//generateListRec(root,potentialVN,vnOut);

}

void Trie::generateListIt(TrieNode* ptr, vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN, vector<unsigned int> &vnOut){

	vector<pair<vector<unsigned int>, vector<unsigned int> > >::iterator it;
	if(ptr == NULL)
		return;
	stack<TrieNode *> stk;
	stk.push(ptr);
	while (!stk.empty()){
		TrieNode *top = stk.top();
		stk.pop();
		vector<TrieNode *> children;
		getChildren(top, children);
		sort(children.begin(), children.end(), sortByLabel());
		for(unsigned int i=0; i<children.size(); i++){
			//cout<<" node "<<top->label<<" child i "<<i<<" label "<<children[i]->label<<" vertexSet ";
			//printV((children[i])->vertexSet);

			stk.push(children[i]);
		}
		//printNodeContent(top);
		sort(top->vertexSet.begin(), top->vertexSet.end());
		int saving = (top->depth - 1) * (top->vertexSet.size() - 1) ;
		if (saving > 0){
			vector<pair<vector<unsigned int>, vector<unsigned int> > > myvect;
			vnOut.push_back(top->label);
			myvect.push_back(make_pair(top->vertexSet, vnOut));
			//printVector(myvect);
			it = search(potentialVN.begin(), potentialVN.end(), myvect.begin(), myvect.end(),CompareList());
			//cout<<"after search \n";
			if(it != potentialVN.end()){
				//if((it->second).size() < saving){
				int fsaving = ((it->second).size() - 1) *((it->first).size() - 1);
				//cout<<"fsaving "<<fsaving<<" saving "<<saving<<" potvn size "<<potentialVN.size()<<" ptr "<<ptr<<endl;
				//cout.flush();	
				
				if(fsaving < saving){
					potentialVN.erase(it);
					//cout<<"erasing verterSet with saving "<<saving<<" fsaving "<<fsaving<<"\n";
				}
			} 
			vnOut.clear();
			vnOut.push_back(top->label);
			//cout<<" vnout "<<top->label<<endl;
			TrieNode *tmp = top->parent;
			while(tmp != NULL){
				//cout<<" parent vnout "<<tmp->label<<endl;
				vnOut.push_back(tmp->label);
				tmp = tmp->parent;	
			}
			sort(vnOut.begin(), vnOut.end());
			potentialVN.push_back(make_pair(top->vertexSet,vnOut)); 
		}
	}
/*
	cout<<"potentialVN ";
	for(it=potentialVN.begin(); it!=potentialVN.end(); it++){
		pair< vector<unsigned int>, vector<unsigned int> > p = *it;
		vector<unsigned int> v1 = p.first;
		for(unsigned int j=0; j<v1.size(); j++){
			cout<<v1[j]<<" ";
		}
		cout<<"- ";
		vector<unsigned int> v2 = p.second;
		for(unsigned int j=0; j<v2.size(); j++){
			cout<<v2[j]<<" ";
		}
		cout<<endl;
	}
*/

}

void Trie::generateListRec(TrieNode* ptr, vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN, vector<unsigned int> &vnOut){
	vector<pair<vector<unsigned int>, vector<unsigned int> > >::iterator it;
	if(vnOut.size() > 17200)
		return;
	if(ptr == NULL)
		return;
	if(ptr != NULL ){
		sort(ptr->vertexSet.begin(), ptr->vertexSet.end());
		int saving = (ptr->depth - 1) * (ptr->vertexSet.size() - 1) ;
		//cout<<" saving "<<saving<<"\n";
		//if(ptr->vertexSet.size() > 0){
		if (saving > 0){
			vector<pair<vector<unsigned int>, vector<unsigned int> > > myvect;
			vnOut.push_back(ptr->label);
			myvect.push_back(make_pair(ptr->vertexSet, vnOut));
			//printVector(myvect);
			
			it = search(potentialVN.begin(), potentialVN.end(), myvect.begin(), myvect.end(),CompareList());
			//cout<<"after search \n";
			//cout.flush();	
			if(it != potentialVN.end()){
				//cout<<"found potential VN erasing\n";
				//cout.flush();	
				//if((it->second).size() < saving){
				int fsaving = ((it->second).size() - 1) *((it->first).size() - 1);
				//cout<<"fsaving "<<fsaving<<" saving "<<saving<<" potvn size "<<potentialVN.size()<<" ptr "<<ptr<<endl;
				//cout.flush();	
				
				if(fsaving < saving){
					potentialVN.erase(it);
					//cout<<"erasing verterSet with saving "<<saving<<" fsaving "<<fsaving<<"\n";
				}
				
				sort(vnOut.begin(), vnOut.end(), sortString());
				//cout<<"after sort vnOut size "<<vnOut.size()<<endl;
				//cout.flush();	
				potentialVN.push_back(make_pair(ptr->vertexSet,vnOut)); 
				//printVector(potentialVN);
				//cout.flush();
			} else {
				//cout<<"adding potential VN \n";
				//cout.flush();	
				sort(vnOut.begin(), vnOut.end(), sortString());
				potentialVN.push_back(make_pair(ptr->vertexSet,vnOut)); 
			}
		}
		generateListRec(ptr->firstChild,potentialVN,vnOut);
		//cout<<"AFTER generateListRec firstchild  parent "<<"\n";
		if(ptr == root)return;
		//cout<<"NEXT parent "<<ptr->parent->label<<"\n";
		vnOut.clear();
		TrieNode *tmp = ptr->parent;
		while(tmp != NULL){
			vnOut.push_back(tmp->label);
			tmp = tmp->parent;	
		}
		sort(vnOut.begin(), vnOut.end());
		generateListRec(ptr->nextSibling,potentialVN,vnOut);
	}
}

void Trie::printVector(vector<pair<vector<unsigned int>, vector<unsigned int> > > & v){
	vector<pair<vector<unsigned int>, vector<unsigned int> > >::iterator it;
	vector<unsigned int>::iterator vit;
	cout<<" PrintVector\n";
	for(it = v.begin(); it != v.end(); it++){
		pair<vector<unsigned int>, vector<unsigned int> > p = *it;
		vector<unsigned int> v1 = p.first;
		vector<unsigned int> v2 = p.second;
		for(vit = v1.begin(); vit != v1.end(); vit++){
			cout<<*vit<<" ";
		}
		cout<<"\n";
		cout<<"vnOut\n";
		for(vit = v2.begin(); vit != v2.end(); vit++){
			cout<<*vit<<" ";
		}
		cout<<"\n";
	}

}

void Trie::printNodeContent(TrieNode *ptr){
	if(ptr == NULL){
		cout<<" NULL\n";
	}
	cout<<ptr->label<<" : ";
	for(vector<unsigned int>::iterator it = ptr->vertexSet.begin(); it != ptr->vertexSet.end(); it++){
		cout<<" "<<*it;
	}
	cout<<" : depth "<<ptr->depth;
	cout<<"\n";
}

void Trie::printTrieRec(TrieNode* ptr){
	if(ptr != NULL){

		if(ptr->parent != NULL)
			cout<<" parent "<<(ptr->parent)->label<<" ";
		cout<<" label "<<ptr->label<<" children "<<ptr->nChildren<<" depth "<<ptr->depth<<" vertexes : ";
		sort(ptr->vertexSet.begin(), ptr->vertexSet.end());
		for(vector<unsigned int>::iterator it = ptr->vertexSet.begin(); it != ptr->vertexSet.end(); it++){
			cout<<" "<<*it;
		}
		cout<<"\n";
		printTrieRec(ptr->firstChild);
		printTrieRec(ptr->nextSibling);
	}
}

void Trie::getChildren(TrieNode *ptr, vector<TrieNode *> &children){
	if(ptr == NULL) return;
	if (ptr->firstChild == NULL) return;
	children.push_back(ptr->firstChild);
	ptr = ptr->firstChild;
	ptr = ptr->nextSibling;
	while(ptr != NULL){
		children.push_back(ptr);
		ptr = ptr->nextSibling;
	}
	return;
}

bool Trie::isEmpty(){
	if(root == NULL)
		return true;
	return false;	
}

void Trie::makeEmpty(){
	if(root != NULL)
		makeEmpty(root);
}

void Trie::makeEmpty(TrieNode * & t){
	if(t != NULL){
		makeEmpty(t->nextSibling);
		makeEmpty(t->firstChild);
		ncount--;
		delete t;
		t = NULL;
		return;
	}
}
	
