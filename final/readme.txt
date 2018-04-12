1. 學號：
B03901001, B03901096
2. 姓名：
廖宜倫, 周晁德
3. 使用之程式語言：< C++ >

4. 使用之編譯器：< GNU g++ >

5. 檔案壓縮方式: zip

6. 各檔案說明：

	readme.txt: 本檔，使用說明及相關資訊
	final_DPnInv.cpp: 程式
	makefile: make file
	 

7. 編譯方式說明： 主程式：final_DPnInv.cpp請在主程式的目錄下，鍵入make指令，即可完成編譯，(makefile原先下的 Optimize 指令為 -O3，若是要改用 -O 請修改同資料匣下的makefile)。在主程式的目錄下會產生一個名為term的執行檔。如果要重新編譯，請先執行 make clean 再執行一次 make。

8. 執行、使用方式說明： 主程式： final_DPnInv.cpp，編譯完成後，在檔案目錄下會產生一個em的執行檔。執行檔的命令格式為:./term <input file name> <output file name>。例如：./term case1 case1.out

9. 執行結果說明（說明執行結果的觀看方法，及解釋各項數據等)：
用8.執行後，會生成case1.out，結果有4行，依序是從北(N)、東(E)、南(S)、西(W)在不同時間點通過十字路口到達不同方向(N, E, S, W)或是沒有車通過(00)

可以用time ./term case1 case1.out測量程式執行的速度。


      

