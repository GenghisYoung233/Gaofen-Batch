const { app, BrowserWindow, ipcMain, dialog, Menu } = require('electron');
const log = require('electron-log');
const path = require('path');

// Configure electron-log
log.transports.file.level = 'info';

function createWindow() {
    const mainWindow = new BrowserWindow({
        width: 800,
        height: 600,
        webPreferences: {
            nodeIntegration: true,
            contextIsolation: false,
        },
        icon: path.join(__dirname, './assets/app_icon.ico'),
        title: 'Gaofen Batch'
    });

    // Set the default menu to null to remove it
    Menu.setApplicationMenu(null);

    mainWindow.loadFile('index.html');
}

app.whenReady().then(createWindow);

app.on('window-all-closed', function () {
    if (process.platform !== 'darwin') app.quit();
});

app.on('activate', function () {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
});

ipcMain.on('open-file-dialog', (event) => {
    dialog.showOpenDialog({
        properties: ['openFile', 'multiSelections'],
        filters: [{ name: 'Tarball Files', extensions: ['tar.gz', 'tar.xz'] }]
    }).then(result => {
        if (!result.canceled) {
            event.sender.send('selected-files', result.filePaths);
        }
    }).catch(err => {
        log.info(err);
    });
});

// New listener for directory selection
ipcMain.on('open-directory-dialog', (event) => {
    dialog.showOpenDialog({
        properties: ['openDirectory'],
        title: '选择输出文件夹'
    }).then(result => {
        if (!result.canceled && result.filePaths[0]) {
            event.sender.send('selected-directory', result.filePaths[0]);
        }
    }).catch(err => {
        log.info(err);
    });
});

ipcMain.on('log', (event, text) => {
    log.info(text);
});
