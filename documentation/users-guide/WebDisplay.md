# Web-based display

In this chapter discussed how web-based display can be created in the ROOT.

## Basics

Idea of web displays is implementation of user interface, which can run remotely
or locally in the web-browsers, fully decoupled from the application code.
For the communication between application and browser websockets are used.


## Creating web-window

**`ROOT::Experimental::TWebWindow`** class is represent window instance, which can be displayed in the browser

```{.cpp}

std::shared_ptr<ROOT::Experimental::TWebWindow> win = ROOT::Experimental::TWebWindowsManager::Instance()->CreateWindow();

// set HTML page which is showed when window displayed
win->SetDefaultPage("file:page.html"); // set

// allow unlimitted user connections to the window (default only 1)
win->SetConnLimit(0);

// configure predefined geometry
win->SetGeometry(300, 300);

```

## Communication

To receive data from the browser, data callback should be assigned.
Callback is invoked when data received from the client or when connection is established.
Normally lambda function is used:

```{.cpp}

win->SetDataCallBack([](unsigned connid, const std::string &msg) {

   printf("Msg:%s from connection:%u\n", msg.c_str(), connid);

});

```

Here **connid** is unique identifier, which assign to each connection when it is established.
There are several predefined messages kinds: **"CONN_READY"** when new connection established and
**"CONN_CLOSED"** when connection is closed by client.
The connection identifier should be used when sending message to the client:

```{.cpp}

// get connection id for the first connection in the list

if (win->NumConnections() > 0) {
   unsigned connid = win->GetConnectionId();
   std::string msg = "Hello, world";
   win->Send(msg, connid);
}

```

## Display window

To display window in the browser, one should call `win->Show()` method.
This will starts new window (or new tab) in the default browser and show content of HTML page,
configured for the window. As argument of `Show()` method one can specify browser kind like
"chromium" or "firefox" or just full path to the program which should be invoked.
With the method `win->GetUrl()` one obtains URL string, which can be typed in the browser address string directly.

Same window can be displayed several times in different browsers or different browser tabs - one only
must allow appropriate number of connections calling ``win->SetConnLimit(3)``

For the local displays **Chromium Embeded Framework (CEF)** is used. It provides functionality
of Chrome web browser in ROOT application without need to create and start real http server.
If CEF was configured correctly, it is enough to call `win->Show("cef")` to display window in CEF.


## Client code

There is no limitations which framework should be used on the client side.
The minimal HTML/JavaScript code, which establish connection with the server, looks like:

``` {.html}
<!DOCTYPE HTML>
<html>
    <head>
        <meta charset="utf-8">
        <title>Web Window Example</title>
        <script src="/jsrootsys/scripts/JSRoot.core.js" type="text/javascript"></script>
    </head>
   <body>
     <div id="main"></div>
     <script>
       JSROOT.connectWebWindow({
          receiver: {
              OnWebsocketOpened: function(handle) {
                  console.log('Connected');
                  handle.Send("Init msg from client");
              },
              OnWebsocketMsg: function(handle, msg) {
                  console.log('Get message ' + msg);
                  document.getElementById("main").innerHTML = msg;
              },
              OnWebsocketClosed: function(handle) {
                 // when connection closed, close panel as well
                 if (window) window.close();
              }
          }
       });
   </script>
   </body>
</html>

```
