#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <wchar.h>

static int show_error(const wchar_t *message)
{
    MessageBoxW(NULL, message, L"Molecule Viewer", MB_OK | MB_ICONERROR);
    return 1;
}

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE previous, wchar_t *arguments, int show)
{
    wchar_t directory[32768];
    wchar_t python[32768];
    wchar_t command[32768];
    const DWORD capacity = ARRAYSIZE(directory);
    DWORD length = GetModuleFileNameW(NULL, directory, capacity);
    if (length == 0 || length == capacity) {
        return show_error(L"Cannot find the application directory.");
    }
    wchar_t *separator = wcsrchr(directory, L'\\');
    if (separator == NULL) {
        return show_error(L"Cannot find the application directory.");
    }
    *separator = L'\0';
    const wchar_t *python_suffix = L"\\runtime\\pythonw.exe";
    const wchar_t *script_suffix = L"\\runtime\\viewer.py";
    if (2 * wcslen(directory) + wcslen(python_suffix) + wcslen(script_suffix) +
        wcslen(arguments) + 7 >= ARRAYSIZE(command)) {
        return show_error(L"The application path or file path is too long.");
    }

    wcscpy(python, directory);
    wcscat(python, python_suffix);
    wcscpy(command, L"\"");
    wcscat(command, python);
    wcscat(command, L"\" \"");
    wcscat(command, directory);
    wcscat(command, script_suffix);
    wcscat(command, L"\" ");
    wcscat(command, arguments);

    STARTUPINFOW startup = {0};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process = {0};
    if (!CreateProcessW(python, command, NULL, NULL, FALSE, 0, NULL, NULL, &startup, &process)) {
        return show_error(L"Extract the full ZIP archive, then run MoleculeViewer.exe again.");
    }

    CloseHandle(process.hThread);
    WaitForSingleObject(process.hProcess, INFINITE);
    DWORD status = 0;
    GetExitCodeProcess(process.hProcess, &status);
    CloseHandle(process.hProcess);
    if (status != 0) {
        return show_error(L"The viewer could not start or exited with an error.");
    }
    return 0;
}
