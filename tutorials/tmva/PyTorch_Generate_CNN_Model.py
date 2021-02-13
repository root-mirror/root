import torch
from torch import nn

# Define model

# Custom Reshape Layer
class Reshape(torch.nn.Module):
    def forward(self, x):
        return x.view(-1,1,16,16)

# CNN Model Definition
net = torch.nn.Sequential(
    Reshape(),
    nn.Conv2d(1, 10, kernel_size=3, padding=1),
    nn.ReLU(),
    nn.BatchNorm2d(10),
    nn.Conv2d(10, 10, kernel_size=3, padding=1),
    nn.ReLU(),
    nn.MaxPool2d(kernel_size=2),
    nn.Flatten(),
    nn.Linear(10*8*8, 256),
    nn.ReLU(),
    nn.Linear(256, 2),
    nn.Sigmoid()
    )

# Construct loss function and Optimizer.
criterion = nn.BCELoss()
optimizer = torch.optim.Adam


def fit(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):
    trainer = optimizer(model.parameters(), lr=0.01)
    schedule, schedulerSteps = scheduler
    best_val = None

    # Setup GPU
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    for epoch in range(num_epochs):
        # Training Loop
        # Set to train mode
        model.train()
        running_train_loss = 0.0
        running_val_loss = 0.0
        for i, (X, y) in enumerate(train_loader):
            trainer.zero_grad()
            X, y = X.to(device), y.to(device)
            output = model(X)
            target = y
            train_loss = criterion(output, target)
            train_loss.backward()
            trainer.step()

            # print train statistics
            running_train_loss += train_loss.item()
            if i % 4 == 3:    # print every 4 mini-batches
                print(f"[{epoch+1}, {i+1}] train loss: {running_train_loss / 4 :.3f}")
                running_train_loss = 0.0

        if schedule:
            schedule(optimizer, epoch, schedulerSteps)

        # Validation Loop
        # Set to eval mode
        model.eval()
        with torch.no_grad():
            for i, (X, y) in enumerate(val_loader):
                X, y = X.to(device), y.to(device)
                output = model(X)
                target = y
                val_loss = criterion(output, target)
                running_val_loss += val_loss.item()

            curr_val = running_val_loss / len(val_loader)
            if save_best:
               if best_val==None:
                   best_val = curr_val
               best_val = save_best(model, curr_val, best_val)

            # print val statistics per epoch
            print(f"[{epoch+1}] val loss: {curr_val :.3f}")
            running_val_loss = 0.0

    print(f"Finished Training on {epoch+1} Epochs!")

    return model


def predict(model, test_X, batch_size=100):
    # Set to eval mode

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    model.eval()


    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    predictions = []
    with torch.no_grad():
        for i, data in enumerate(test_loader):
            X = data[0].to(device)
            outputs = model(X)
            predictions.append(outputs)
        preds = torch.cat(predictions)

    return preds.cpu().numpy()


load_model_custom_objects = {"optimizer": optimizer, "criterion": criterion, "train_func": fit, "predict_func": predict}

# Store model to file
m = torch.jit.script(net)
torch.jit.save(m,"PyTorchModelCNN.pt")
