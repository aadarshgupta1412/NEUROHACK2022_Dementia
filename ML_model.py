#!/usr/bin/env python3

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import random
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader
import torchvision.datasets as datasets
import torchvision.transforms as transforms
import sys

#load data

path1=os.getcwd()
path0 = os.path.dirname(path1)
path2=path0+'/out'

df=pd.read_stata(path0+'/H_DAD_w1a3.dta')

df.set_index('prim_key', inplace=True)

ft_lst = ['r1borient',
'r1bexefu',
'r1blangf',
'r1bmemory',
'r1bvsp',
'r1nmemimm',
'r1nmemdel',
'r1nmemrec',
'r1nreason',
'r1natnspd',
'r1sgcp',
'r1hmse_scorz',
'r1word_totaz',
'r1word_dz',
'r1wre_scorez',
'r1log_recoz',
'r1bm_immexz',
'r1bm_reclexz',
'r1verbalz',
'r1csid_scorz',
'r1rv_scorez',
'r1cog_totalz']

out_lst = ['r1inf_strok',
'r1inf_parkn',
'r1inf_alzhe',
'r1inf_memry'
]




##create fully connected network

class NN(nn.Module): #
    def __init__(self, input_size, num_classes):
        super(NN, self).__init__()
        self.fc1 = nn.Linear(input_size, 50)
        self.fc2 = nn.Linear(50, num_classes)
        
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x


#model = NN(784 ,10)
#x = torch.randn(64, 784)
#print(model(x).shape)
#sys.exit()

#set device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#hpyerparameters
input_size = 22
num_classes = 16
learning_rate = 0.001
batch_size = 64
num_epochs = 1



##Load data
df_features = df[ft_lst]
df_out=df[out_lst]

mask=np.ones(df_features.shape[0],dtype=bool)
mask[0:round(df_features.shape[0]*0.2)]=False
random.shuffle(mask)




train_dataset = torch.tensor(df[mask].values.astype(np.float32))
train_loader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True)
test_dataset = torch.tensor(df[~mask].values.astype(np.float32))
test_loader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=True)



#Initialize network

model = NN(input_size = input_size, num_classes=num_classes).to(device)



#Loss and optimizer

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=learning_rate)


#Train Network

for epcoh in range(num_epochs):
    for batch_idx, (data, targets) in enumerate(train_loader):
        
        #Get data to cuda if possible
        data = data.to(device=device)
        targets = targets.to(device=device)
        
        #get to correct shape
        data = data.reshape(data.shape[0], -1)
        
        
        #forward
        scores = model(data)
        loss = criterion(scores, targets)
        
        # backward 
        optimizer.zero_grad()
        loss.backward()
        
        #gradieant descent or adam step
        optimizer.step()
        



##Check accuracy on training & test to see how good our model is

def check_accuracy(loader, model):
    
    num_correct = 0 
    num_samples = 0
    model.eval()
    
    with torch.no_grad():
        for x,y in loader:
            x= x.to(device=device)
            y=y.to(device=device)
            x=x.reshape(x.shape[0],-1)
            
            scores = model(x)
            _, predictions = scores.max(1)
            num_correct += (predictions == y).sum()
            num_samples += predictions.size(0)
            
        print(f'Got {num_correct}/{num_samples} with accuracy {float(num_correct)/float(num_samples)*100:.2f}]')
    model.train()
    #return acc

check_accuracy((train_loader), model)
check_accuracy(test_loader, model)