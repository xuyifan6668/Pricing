@echo off
rem Developer shell helper for Windows: sets up MSVC, CMake, Ninja and vcpkg, then runs the given
rem command in the repo root. Examples:
rem   tools\dev.cmd cmake --preset msvc-debug
rem   tools\dev.cmd cmake --build --preset msvc-debug
rem   tools\dev.cmd ctest --preset msvc-debug
rem   tools\dev.cmd cmd            (interactive shell with everything on PATH)
setlocal
if "%VCPKG_ROOT%"=="" set "VCPKG_ROOT=C:\dev\vcpkg"
set "PATH=C:\Program Files\CMake\bin;%LOCALAPPDATA%\Microsoft\WinGet\Packages\Ninja-build.Ninja_Microsoft.Winget.Source_8wekyb3d8bbwe;%LOCALAPPDATA%\Microsoft\WinGet\Links;%PATH%"
set "VSDEV=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VSDEV%" set "VSDEV=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
call "%VSDEV%" >nul 2>&1
cd /d "%~dp0.."
%*
endlocal & exit /b %ERRORLEVEL%
