import torch
from torch import nn

# Define model
model = nn.Sequential(
                nn.Linear(4, 64),
                nn.ReLU(),
                nn.Linear(64, 2),
                nn.Softmax(dim=1))

# Construct loss function and Optimizer.
criterion = torch.nn.MSELoss()
optimizer = torch.optim.SGD


def fit(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):
    trainer = optimizer(model.parameters(), lr=0.01)
    schedule, schedulerSteps = scheduler
    best_val = None

    for epoch in range(num_epochs):
        # Training Loop
        # Set to train mode
        model.train()
        running_train_loss = 0.0
        running_val_loss = 0.0
        for i, (X, y) in enumerate(train_loader):
            trainer.zero_grad()
            output = model(X)
            train_loss = criterion(output, y)
            train_loss.backward()
            trainer.step()

            # print train statistics
            running_train_loss += train_loss.item()
            if i % 32 == 31:    # print every 32 mini-batches
                print(f"[{epoch+1}, {i+1}] train loss: {running_train_loss / 32 :.3f}")
                running_train_loss = 0.0

        if schedule:
            schedule(optimizer, epoch, schedulerSteps)

        # Validation Loop
        # Set to eval mode
        model.eval()
        with torch.no_grad():
            for i, (X, y) in enumerate(val_loader):
                output = model(X)
                val_loss = criterion(output, y)
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


def predict(model, test_X, batch_size=32):
    # Set to eval mode
    model.eval()
   
    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    predictions = []
    with torch.no_grad():
        for i, data in enumerate(test_loader):
            X = data[0]
            outputs = model(X)
            predictions.append(outputs)
        preds = torch.cat(predictions)
   
    return preds.numpy()


load_model_custom_objects = {"optimizer": optimizer, "criterion": criterion, "train_func": fit, "predict_func": predict}

# Store model to file
m = torch.jit.script(model)
torch.jit.save(m,"PyTorchModelClassification.pt")

