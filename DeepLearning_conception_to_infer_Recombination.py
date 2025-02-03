

import pandas as pd
path='/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1'
#data = pd.read_csv(path+'/Candiate_recombinations/cp22_IBDclean_w80_cov.Recomb')
data = pd.read_csv(path+'/Candiate_recombinations/cp22_IBDclean_w80_covariable_1CPbyRecomb.Recomb')

#Filter NA
data22 = data.dropna(subset=['index.G', 'index.D'])

#keep CPs >0.5
data22 = data22[data22['PROB'] >= 0.5]

import random
# Fix the seed for reproducibility
random.seed(42)

# Sample 100 items from the list
#data22 =  data22.iloc[random.sample(range(data22.shape[0]), 1_000), :] #<---------------------

# Model.....
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import torch.nn.init as init
import pandas as pd
import numpy as np


class WaterDataset(Dataset):
    def __init__(self, dfpd):
        super().__init__()
        # Load data to pandas DataFrame
        #df = pd.read_csv(csv_path)
        # Convert data to a NumPy array and assign to self.data
        self.data = dfpd

    # Implement __len__ to return the number of data samples
    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, idx):
        features = self.data[idx, :-1]
        # Assign last data column to label
        label = self.data[idx, -1]
        return features, label


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        # Define the three linear layers
        self.fc1 = nn.Linear(15,150)
        self.fc2 = nn.Linear(150,100)
        self.fc3 = nn.Linear(100,50)
        self.fc4 = nn.Linear(50,1)

        # Add two batch normalization layers
        self.bn1 =nn.BatchNorm1d(150)
        self.bn2 =nn.BatchNorm1d(100)
        self.bn3 =nn.BatchNorm1d(50)


        # Apply He initialization
        init.kaiming_uniform_(self.fc1.weight)
        init.kaiming_uniform_(self.fc2.weight)
        init.kaiming_uniform_(self.fc3.weight)
        init.kaiming_uniform_(self.fc4.weight, nonlinearity="sigmoid")


    def forward(self, x):
        # Pass x through linear layers adding activations
        x=self.fc1(x)
        x = nn.functional.elu(self.bn1(x))
        x=self.fc2(x)
        x = nn.functional.elu(self.bn2(x))
        x=self.fc3(x)
        x = nn.functional.elu(self.bn3(x))
        x = nn.functional.sigmoid(self.fc4(x))
        return x


# Convert the NumPy array to a PyTorch tensor
datachr22 = torch.tensor(data22.iloc[:,2:].values)
dataset_train = WaterDataset(datachr22)

# Create a DataLoader based on dataset_train
dataloader_train = DataLoader(dataset_train,  batch_size=(4*64),  shuffle=True)

# Get a batch of features and labels
features, labels = next(iter(dataloader_train))
#print(features, labels)

# Get unique labels and their counts
unique_labels, counts = torch.unique(labels, return_counts=True)

# Combine them into a dictionary
table = dict(zip(unique_labels.tolist(), counts.tolist()))
print("unique labels and their counts for the last batch")
print(table)


# Function to save the model _1CPbyRecomb
#filename=path+'/TunningRDL/best_model.pth'
filename=path+'/TunningRDL/best_model_1CPbyRecomb.pth'
def save_checkpoint(model, optimizer, epoch, loss, filename=filename):
    checkpoint = {
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'epoch': epoch,
        'loss': loss
    }
    torch.save(checkpoint, filename)



# Function to load the model
def load_checkpoint( filename=filename):
    checkpoint = torch.load(filename)
    model=(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
    learning=(checkpoint['optimizer_state_dict'])
    epoch = checkpoint['epoch']
    loss = checkpoint['loss']
    return model, optimizer,learning, epoch, loss,

def train_model(optimizer, net, num_epochs):
    best_loss = float('inf')
    best_epoch = -1
    best_model_weights = None

    for epoch in range(num_epochs):
        for features, labels in dataloader_train:
            optimizer.zero_grad()

            # Convert input data to the same dtype as the weight parameters
            features = features.to(torch.float)  # Adjust the data type as needed

            outputs = net(features)
            labels = labels.to(torch.float)

            loss = criterion(outputs, labels.view(-1, 1))

            loss.backward()
            optimizer.step()

        # Save the model weights at each epoch
        #torch.save(net.state_dict(), f'model_epoch_{epoch}.pth')

        # Check if this is the best model so far
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_epoch = epoch
            best_model_weights = net.state_dict()
            save_checkpoint(net, optimizer, best_epoch, best_loss)

        if epoch % 100 == 0:
          print('Epoch {}: Loss = {}'.format(epoch, loss.item()))


          if loss.item() > best_loss:
            #print("loss.item() > best_loss")
            #Load the best model
            model, optimizer,learning,  epoch, loss = load_checkpoint()
            #model, _ = load_checkpoint('best_model.pth')
            learning_rate=learning['param_groups'][0]['lr']
            net.load_state_dict(model)
            optimizer = optim.Adam(net.parameters(), lr=learning_rate) #/1.2
          else:
            #print("loss.item() <= best_loss")
            learning_rate=optimizer.param_groups[0]['lr']
            optimizer = optim.Adam(net.parameters(), lr=learning_rate) #*1.1


    # Save the best model weights
    #torch.save(best_model_weights, 'best_model3.pth')
    print('Best Epoch: {}, Best Loss: {}'.format(best_epoch, best_loss))

    return best_epoch, best_loss



import torch
import torch.optim as optim
import torch.nn as nn

criterion = nn.BCELoss()
net = Net()
optimizer = optim.Adam(net.parameters())

import os
print("Weights before loading:")
print(list(net.parameters())[0][0][:5])
if os.path.isfile(filename):
    print("Already have one best model")
    #Load the saved weights
    model, optimizer, learning, epoch, loss = load_checkpoint(filename)
    net.load_state_dict(model)
    learning_rate = 1e-6 #learning['param_groups'][0]['lr']

else:
    print("No best model saved, starting from zero ...")
    optimizer = optim.Adam(net.parameters(), lr=1e-6)
    learning_rate = optimizer.param_groups[0]['lr']


print("Weights after loading:")
print(list(net.parameters())[0][0][:5])

print(f"Learning Rate befor: {learning_rate}")

print("Training model ....")
train_model(optimizer, net, num_epochs=5000) #<----------------------------------

#for param_group in optimizer.param_groups:
#    learning_rate = param_group['lr']
learning_rate = optimizer.param_groups[0]['lr']
print(f"Learning Rate After this train : {learning_rate}")

from torchmetrics import Accuracy, ConfusionMatrix
# Set up binary accuracy metric
acc = Accuracy(task="binary")
confmat = ConfusionMatrix(task="binary", num_classes=2)

all_preds = []
all_labels = []

net.eval()
with torch.no_grad():
    for features, labels in dataloader_train:  # or dataloader_test
        # Get predicted probabilities for test data batch
        features = features.to(torch.float32)
        outputs = net(features)
        labels = labels.to(torch.float32)
        preds = (outputs > 0.5).float()

        # Accumulate predictions and labels
        all_preds.append(preds)
        all_labels.append(labels.view(-1, 1))

# Concatenate all predictions and labels
all_preds = torch.cat(all_preds)
all_labels = torch.cat(all_labels)

# Compute accuracy and confusion matrix
accuracy = acc(all_preds, all_labels)
confusion_matrix = confmat(all_preds, all_labels)

print(f"Accuracy: {accuracy}")
print(f"Confusion Matrix:\n{confusion_matrix}")

