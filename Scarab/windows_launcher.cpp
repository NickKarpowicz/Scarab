#define WIN32_LEAN_AND_MEAN
#include <windows.h>

int wmain(void)
{
    wchar_t exe_path[MAX_PATH];
    GetModuleFileNameW(NULL, exe_path, MAX_PATH);
    PathRemoveFileSpecW(exe_path);

    wchar_t bin_path[MAX_PATH];
    swprintf(bin_path, L"%s\\bin", exe_path);

    SetEnvironmentVariableW(L"GTK_EXE_PREFIX", bin_path);
    SetEnvironmentVariableW(L"GTK_DATA_PREFIX", exe_path);
    SetEnvironmentVariableW(L"GSETTINGS_SCHEMA_DIR",
        L"%s\\share\\glib-2.0\\schemas", exe_path);

    SetDllDirectoryW(L"%s\\lib");

    wchar_t command_line_call[MAX_PATH * 2];
    swprintf(command_line_call, L"\"%s\\myapp.exe\"", bin_path);
    STARTUPINFOW si = { sizeof(si) };
    PROCESS_INFORMATION pi;
    CreateProcessW(NULL, command_line_call, NULL, NULL, FALSE, 0,
                   NULL, NULL, &si, &pi);
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD rc;
    GetExitCodeProcess(pi.hProcess, &rc);
    return (int)rc;
}
