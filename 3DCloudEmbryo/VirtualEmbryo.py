#!/usr/bin/env python
# -*- coding: latin-1 -*-

import sys
import os
import time
import bz2
import httplib,urllib
import hashlib
import json

def addslashes(s):
    d = {'"':'\\"', "'":"\\'", "\0":"\\\0", "\\":"\\\\"}
    return ''.join(d.get(c, c) for c in s)


def tryParseInt(value):
    try:
        return int(value), True
    except ValueError:
        return value, False

class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

ss=" --> "
def strblue(strs):
    return bcolors.BLUE+strs+bcolors.ENDC
def strred(strs):
    return bcolors.RED+strs+bcolors.ENDC
def strgreen(strs):
    return bcolors.BOLD+strs+bcolors.ENDC


class VirtualEmbryo:

    #COMMUNICATION PARAMETERS TO 3D CLOUD EMBRYOS
    url="embryon3d.crbm.cnrs.fr"
    port=80
    Hash="ascidie12:"
    headers = {"Content-type": "application/x-www-form-urlencoded","Accept": "text/plain"}
    
    def __init__(self,login,passwd): # CONSTRUCTATORRRR
        
        self.id_people=-1
        self.id_scenario=-1
        self.id_scenario_owner=-1
        self.id_dataset=-1
        self.id_dataset_owner=-1
        self.id_selection=-1
        self.id_selection_owner=-1
        self.minTime=-1
        self.maxTime=-1
        self.timepoint=[]
        self.lineage="";
        self.login=login
        self.passwd=passwd
        self.selection="{"
        self.connect()

        
    #CONNECTION
    def connect(self):
        conn = httplib.HTTPConnection(self.url,80)
        params = urllib.urlencode({'hash': self.Hash, 'login': self.login, 'passwd': hashlib.md5(self.passwd).hexdigest()})
        conn.request("POST", "/api.php",params, self.headers)
        response=conn.getresponse()
        if response.status==200:
            data=response.read()
            if len(data)>4 and data[0:5]=="ERROR":
                print data
                return False
            else:
                self.id_people = data
                print strblue(self.login+' connected to 3D Cloud Embryos with id '+self.id_people)
                return True
        else:
            print strred('CONNECTION ERROR '+str(response.status)+" "+response.reason)
        return False

    def isConnected(self):
        if self.id_people==-1:
            print strred(' ERROR : You are not connected ')
            quit()
            return False
        return True
    def request(self,param):
        if self.isConnected():
            conn = httplib.HTTPConnection(self.url,80)
            #Add Default Parameters
            param['hash']=self.Hash
            param['id_people']=self.id_people
            params = urllib.urlencode(param)
            conn.request("POST", "/api.php",params, self.headers)
            response=conn.getresponse()
            if response.status==200:
                return response.read()
            else:
                print strred('CONNECTION ERROR '+str(response.status)+" "+response.reason)
                quit()

    #LOOK FOR SOMEONE IN THE DATABASE
    def getGuyByID(self,id_guy): # RETURN NAME + SURNAME + LOGIN FORM SPECIFIC ID
        data=self.request({"who":id_guy})
        if data=="[]":
            print strblue(ss+str(id_guy)+" is unkown")
            return -1
        else:
            #print str(id_guy) +" is "+data
            return id_guy
    def getGuyByName(self,name): # RETURN NAME + SURNAME + LOGIN FORM SPECIFIC ID
        data=self.request({"whoname":name})
        if data=="[]":
            print strblue(ss+str(name)+" is unkown")
            return -1
        else:
            dataset=json.loads(data)[0]
            #print str(name) +" id is "+str(dataset[0])
            return int(dataset[0])
    
    #SHARE
    def share(self,base,id_base,who): #SHARE WHATEVER
        if self.isConnected():
            id_who=-1
            if who=="all":
                id_who=0
            else:
                id_who=self.getGuyByName(who)
                if id_who==-1:
                    return False
                else:
                    data=self.request({"share":id_who,"base":base,"id_base":id_base})
                    print  ss+base +' '+str(id_base)+ " share with " +who + " on id "+data
    def deleteShare(self,base,id_base,who): #DELETE A SPECIFIC SHARE
        if self.isConnected():
            id_who=-1
            if who=="all":
                id_who=0
            else:
                id_who=self.getGuyByName(who)
                if id_who==-1:
                    return False
                else:
                    data=self.request({"deleteshare":id_who,"base":base,"id_base":id_base})
                    if data=="":
                        print  ss+"remove shared " +base +' '+str(id_base)+ " with " +who + " on id "+data
                    else:
                        print strred(' ERROR : cannot delete this share '+data)
    def deleteShareAll(self,base,id_base): #DELETE A ALL SHARE
        if self.isConnected():
            data=self.request({"deleteshareall":0,"base":base,"id_base":id_base})
            if data=="":
                print  ss+"remove all shared " +base +' '+str(id_base)
            else:
                print strred(' ERROR : cannot delete this share '+data)
    
    
    #DATASET
    def isDataSet(self):
        if not self.isConnected():
            return False
        if self.id_dataset==-1:
            print strgreen(ss+'you first have to select a dataset')
            return False
        return True
    def ownDataSet(self):
        if not self.isDataSet():
            return False
        if self.id_dataset_owner!=self.id_people:
            print strgreen(ss+'you are not the owner of this dataset, ask '+self.getGuy(self.id_dataset_owner))
            return False
        return True  
    def initTimePoint(self,minTime,maxTime): #INTERNAL FUNCTION TO INITIALISE TIME POINT
        self.minTime=int(minTime)
        self.maxTime=int(maxTime)
        self.timepoint = ["" for x in range(self.maxTime+1)]
    def createDataSet(self,dataname,minTime,maxTime): #CREATE A NEW DATA SET
        if self.isConnected():
            data=self.request({"createdataset":dataname,"minTime":minTime,"maxTime":maxTime})
            self.id_dataset_owner=self.id_people
            ids,ok=tryParseInt(data)
            if not ok: 
                print strred(' ERROR : Dataset not created '+data)
            else :
                self.id_dataset=ids
                self.id_dataset_owner=self.id_people
                self.initTimePoint(minTime,maxTime)
                print ss+"your id dataset is "+str(self.id_dataset)
    def parseDataSet(self,data): #Parse id,name,minTime,maxTime,id_people to dataset structure
        if data=="[]" or data=="":
            print strred(ss+'dataset not found')
        else:
            dataset=json.loads(data)[0]
            ids,ok=tryParseInt(dataset[0])
            if not ok: 
                print strgreen(ss+'dataset not found '+data)
            else: 
                name=dataset[1]
                self.initTimePoint(dataset[2],dataset[3])
                self.id_dataset_owner=dataset[4]
                self.id_dataset=ids
                print ss+'found dataset '+name+' with id ' + str(self.id_dataset)+' from '+str(self.minTime)+' to ' +str(self.maxTime)+' owned by '+str(self.id_dataset_owner)
    def selectDataSetById(self,ids): #SELECT A DATASET BY ID
        if self.isConnected():
            self.id_dataset=-1
            data=self.request({"dataset":ids})
            if data=="[]" or data=="": #Data set not find ->LOOK IN SHARE DATA SET
                data=self.request({"sharedataset":ids})
            self.parseDataSet(data)           
    def selectDataSetByName(self,name): #SELECT A DATASET BY NAME
        if self.isConnected():
            self.id_dataset=-1
            data=self.request({"datasetname":name})
            if data=="[]" or data=="": #Data set not find ->LOOK IN SHARE DATA SET
                data=self.request({"sharedatasetname":name})
            self.parseDataSet(data)
    def shareDataSet(self,who): #SHARE A DATA SET WITH A GUY GIVEN NAME + SURNAME
        if self.ownDataSet():
            self.share("dataset",self.id_dataset,who)
    def unShareDataSetWith(self,who):#DELETE A SHARE A DATA SET WITH A GUY GIVEN NAME + SURNAME
        if self.ownDataSet():
            self.deleteShare("dataset",self.id_dataset,who)
    def unShareDataSetWithAll(self):#DELETE A SHARE A DATA SET WITH A GUY GIVEN NAME + SURNAME
        if self.ownDataSet():
            self.deleteShareAll("dataset",self.id_dataset)
    def deleteDataSet(self): #COMPLETE DELETE OF A DATASET
        if self.ownDataSet():
            data=self.request({"deletedataset":self.id_dataset})
            if data=="": 
                print ss+'dataset deleted'
                self.id_dataset=-1
                self.id_dataset_owner=-1
                self.minTime=-1
                self.maxTime=-1
            else:
                print strred(' ERROR : '+data)
    def clearDataSet(self): # CLEAR ALL TIME POINT 
        if self.ownDataSet():
            data=self.request({"cleardataset":self.id_dataset})
            if data=="": 
                print ss+'dataset cleared'
            else:
                print strred(' ERROR : '+data)
    def addCell(self,t,id_cell,id_mother,obj):#ADD A CELL CELL
        if self.ownDataSet():
            obj="g cell_"+str(id_cell)+"\n"+obj; #ADD THE GGROUP
            if id_mother!=-1:
                self.lineage+=str(t-1)+";"+str(id_mother)+';'+str(id_cell)+"\n";
            
            self.timepoint[t]+=obj+"\n"; 
    def uploadTimePoint(self,t): #UPLOAD TIME POINT IN DATASET
        if self.ownDataSet():
            data=self.request({"uploadtimepoint":self.id_dataset,"t":t,"data":bz2.compress(self.timepoint[t])})
            ids,ok=tryParseInt(data)
            if not ok: 
                print strred(' ERROR : time point not uploaded '+data)
            else :
                print ss+"time point "+str(t)+" upload with id "+str(ids)
    def uploadInfos(self,infos,field,type): #UPLOAD A INFORMATION ON THE DATASET (coulbe be lineage or name etc..)
        if self.ownDataSet():
            data=self.request({"uploadinfos":self.id_dataset,"infos":infos,"type":type,"field":bz2.compress(field)})
            ids,ok=tryParseInt(data)
            if not ok: 
                print strred(' ERROR : '+infos+' not uploaded '+data)
            else :
                print ss+infos+" upload with id "+str(ids)
    def deleteInfos(self,infos): #DELETE AN UPLOADED INFOS
        if self.ownDataSet():
            data=self.request({"deleteinfos":self.id_dataset,"infos":infos})
            ids,ok=tryParseInt(data)
            if not ok: 
                print strred(' ERROR : '+infos+' '+data)
            else :
                print ss+infos+" with id "+str(ids)+" deleted"
    def uploadLineage(self): #UPLOAD THE LINEAGE FILE (AS INFOS)
        if self.ownDataSet():
            self.uploadInfos("Lineage",self.lineage,'string')
    def deleteLineage(self):
        if self.ownDataSet():
            self.deleteInfos("Lineage")
    def uploadDataset(self): #UPLOAD ALL TIME POINTS AND THE ASSOCIATED LINEAGE 
        if self.ownDataSet():
            for t in range(self.minTime,self.maxTime+1):
                self.uploadTimePoint(t)
            self.uploadLineage()



    '''    
    #SCENARIO
    def isScenario(self):
        if not self.isConnected():
            return False
        if self.id_scenario==-1:
            print ' You first have to select a scenario'
            return False
        return True
    def ownScenario(self):
        if not self.isScenario():
            return False
        if self.id_scenario_owner!=self.id_people:
            print ' You are not the owner of this scenario, ask '+self.getGuy(self.id_scenario_owner)
            return False
        return True
    def createScenario(self,scenarname): #CREATE A NEW SCENARIO
        if self.isDataSet():
            try:
                self.cur.execute('INSERT  INTO scenario (id_people,date,name,id_dataset) VALUES ('+str(self.id_people)+',"'+time.strftime("%Y-%m-%d %H:%S:%I")+'","'+scenarname+'",'+str(self.id_dataset)+')');
                self.db.commit()
                self.id_scenario=self.cur.lastrowid
                self.id_scenario_owner=self.id_people
                print "--> Your id scenario is "+str(self.id_scenario)
            except MySQLdb.Error, e:
                self.catch(e)
    def selectScenarioById(self,ids): #SELECT A SCENARIO BY ID
        if self.isConnected():
            try:
                self.cur.execute('SELECT id,id_dataset FROM scenario WHERE id='+str(ids)+' and id_people='+str(self.id_people));
                for row in self.cur.fetchall():
                    self.id_scenario=row[0]  
                    self.id_dataset=row[1]  
                    self.id_scenario_owner=self.id_people
                if self.id_scenario==-1: #LOOK IN SHARE SCENARIO 
                    self.cur.execute('SELECT id_people FROM share WHERE base="scenario" and id_base='+str(ids)+' and ( id_who='+str(self.id_people)+ ' or id_who=0)' );
                    for row in self.cur.fetchall():
                        temp_id_scenario_owner=row[0]
                        self.cur.execute('SELECT id,id_dataset FROM dataset WHERE id='+str(ids)+' and id_people='+str(temp_id_scenario_owner));
                        for rows in self.cur.fetchall():
                            self.id_scenario=rows[0]  
                            self.id_dataset=rows[1] 
                            self.id_scenario_owner=temp_id_scenario_owner 
                    if self.id_scenario!=-1:
                         print "--> Your shared id scenario is "+str(self.id_scenario)
                    else :
                        print ' Scenario not found' 
                else:
                    print "--> Your id scenario is "+str(self.id_scenario)
            except MySQLdb.Error, e:
                self.catch(e)
    def selectScenarioByName(self,name): #SELECT A SCENARIO BY NAME
        if self.isConnected():
            try:
                self.cur.execute('SELECT id,id_dataset FROM scenario WHERE LOWER(name) like "%'+name.lower()+'%" and id_people='+str(self.id_people));
                for row in self.cur.fetchall():
                    if self.id_scenario==-1:
                        self.id_scenario=row[0] 
                        self.id_dataset=row[1]  
                if self.id_scenario==-1: 
                    print ' Scenario not found ' 
                else:
                    print "--> Your id scenario is "+str(self.id_scenario)
            except MySQLdb.Error, e:
                self.catch(e)
    def deleteScenario(self): #DELETE ALL THE SCENARIO FOR THIS USER
        if self.ownScenario():
            try:
                self.cur.execute('DELETE FROM scenario WHERE id='+str(self.id_scenario)+' and id_people='+str(self.id_people));
                self.cur.execute('DELETE FROM action WHERE id_scenario='+str(self.id_scenario));
                self.deleteShareAll("scenario",self.id_scenario);
                self.db.commit()
                self.id_scenario=-1 
            except MySQLdb.Error, e:
                self.catch(e)
    def shareScenario(self,who):
        if self.ownScenario():
            self.share("scenario",self.id_scenario,who)
    def playScenario(self):
        print 'not yet implemented'

    #FRAME IN SCENARIO
    def deleteFrames(self): #Delete all actions from a scenario
        if self.ownScenario():
            try:
                self.cur.execute('DELETE FROM  action WHERE id_scenario='+str(self.id_scenario));
                self.db.commit()
                print ' --> delete all actions for scenario '+str(self.id_scenario)
            except MySQLdb.Error, e:
                self.catch(e)
    def addFrame(self,time,position,rotation,scale,framerate,interpolate): #Add a action to the scenario
        if self.ownScenario():
            try:
                #Rotation format : (-0.2, -0.2, -0.1, 0.9)
                #Position format : (0.0, 0.0, 0.0)
                try:
                    positionS='('+str(position[0])+','+str(position[1])+','+str(position[2])+')';
                    rotationS='('+str(rotation[0])+','+str(rotation[1])+','+str(rotation[2])+','+str(rotation[3])+')';
                except IndexError:
                    print "Position Format = array [ x , y , z ]"
                    print "Rotation Format = array [ x , y , z , w]"
                id_last_action=0
                self.cur.execute('SELECT idx FROM action WHERE id_scenario='+str(self.id_scenario)+' order by idx');
                for row in self.cur.fetchall():
                    id_last_action=int(row[0])
                id_last_action+=1
                self.cur.execute('INSERT INTO action (id_scenario,idx,time,position,rotation,scale,framerate,interpolate) VALUES ('+str(self.id_scenario)+','+str(id_last_action)+','+str(time)+',"'+positionS+'","'+rotationS+'","'+str(scale)+'",'+str(framerate)+','+str(interpolate).lower()+')');
                self.db.commit()
                print ' --> add action '+str(id_last_action)
            except MySQLdb.Error, e:
                self.catch(e)

    '''


    #SELECTION
    def isSelection(self):
        if not self.isDataSet():
            return False
        if self.id_selection==-1:
             print strgreen(ss+'you first have to create or select a selection')
             return False
        return True
    def ownSelection(self):
        if not self.isSelection():
            return False
        if self.id_selection_owner!=self.id_people:
            print strgreen(ss+'you are not the owner of this selection, ask '+self.getGuy(self.id_selection_owner))
            return False
        return True
    def createSelection(self,name): #CREATE A NEW SELECTION
        if self.isDataSet():
            data=self.request({"createselection":self.id_dataset,"name":name})
            ids,ok=tryParseInt(data)
            if not ok: 
                print strred(' ERROR : cannot create selection '+data)
            else :
                self.id_selection=ids
                self.id_selection_owner=self.id_people
                print ss+"selection "+name+" created with id "+str(self.id_selection)
    def getSelectionById(self,ids): #GET A SELECTION NUMBER AND RETURN LIST OF CELLS
        cells=""
        if self.isDataSet():
            data=self.request({"getselectionbyid":self.id_dataset,"id":ids})
            print data
            #if data=="[]":
            #    data=self.request({"getsharedselectionbyid":self.id_dataset,"id":ids})
        return cells
    def getSelectionByName(self,name): #GET A SELECTION NAME AND RETURN LIST OF CELLS
        cells=""
        if self.isDataSet():
            data=self.request({"getselectionbyname":self.id_dataset,"name":name})
            print data
        return cells
    def addSelectionToCell(self,id_cell,selection_number):#ADD A CELL WITH SELECTION NUMBER
        if self.ownSelection():
            self.selection+='"'+str(id_cell)+'":"'+str(selection_number)+'",'
    def updateSelection(self):
        if self.ownSelection():
            if self.selection[0]!="{":
                self.selection="{"+self.selection
            if self.selection[len(self.selection)-1]==",":
                self.selection=self.selection[0:len(self.selection)-1]
            self.selection+="}"
            data=self.request({"uploadselection":self.id_selection,"cells":self.selection})
            self.selection="{"
            if data!='':
                print strred(' ERROR : cannot upload this selection '+data)
            else :
                print ss+" selection "+str(self.id_selection)+" uploaded"

    def deleteSelection(self):
        if self.ownSelection():
            data=self.request({"deleteselection":self.id_selection})
            if data!='':
                print strred(' ERROR : cannot delete this selection '+data)
            else :
                print ss+" selection "+str(self.id_selection)+" deleted"
                self.id_selection=-1
                self.id_selection_owner=-1
    def deleteAllSelections(self):
        if self.ownDataSet():
            data=self.request({"deleteallselections":self.id_selection})
            if data!='':
                print strred(' ERROR : cannot delete all selections '+data)
            else :
                print ss+" all selection deleted for this dataset "+self.id_dataset
                self.id_selection=-1
                self.id_selection_owner=-1
    def shareSelection(self,who):
        if self.ownSelection():
            self.share("selection",self.id_selection,who)
    def deleteShareSelection(self,who):
        if self.ownSelection():
            self.deleteShare("selection",self.id_selection,who)
    

    #UNITY
    def play(self,action,value):#GENERIC FUNCTION TO PLAY A ACTION IN UTNITY
        data=self.request({"play":action,"value":value})
        ids,ok=tryParseInt(data)
        if not ok: 
            print strred(' ERROR : cannot play '+action+' '+str(value)+' -> '+data)

    def playDataSet(self): #LOAD DATASET IN UNITY
        self.play("dataset",self.id_dataset)  
    def playTime(self,id_time): #MOVE TO A SPECIFIC TIME POINT
        self.play("time",id_time)
    def playRotate(self,rotate): #ROTATE THE EMBRYO WITH  QUATERNION FORMAT (x,y,z,w) 
        self.play("rotate",rotate)
    def playMove(self,move): #MOVE THE EMBRYO WITH  VECTOR3 FORMAT (x,y,z) 
        self.play("move",move)
    def playSelection(self): #LOAD A SELECTION 
        if self.isSelection():
            self.play("selection",self.id_selection)
    def playCancelSelection(self): #RESET ALL SELECTION 
        self.play("resetselection","-1")
    def playCancelSelectionTime(self,id_time):  #RESET ALL SELECTION AT A SPECIFIC TIME 
        self.play("resetselection",id_time)
    def playHideSelection(self,id_hide):  #HIDE CELLS WITH A SPECIFC NUMBER INSIDE THE SLECTION  
        self.play("hideselection",id_hide)
    def playShowSelection(self,id_show):  #SHOW CELLS WITH A SPECIFC NUMBER INSIDE THE SLECTION  
        self.play("showselection",id_show)


   

