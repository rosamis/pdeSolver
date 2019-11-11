#!/bin/bash
# #remove todos os arquivos gerados anteriormente para o novo grafico
export PATH="/home/soft/likwid/bin:/home/soft/likwid/sbin:${PATH}"
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

echo "MODO: PERFORMANCE"

opcao=5
while [ $opcao -ne 7 ]
do
	echo " ==================================="
    echo "|                                   |"
	echo "| Escolha uma opção:                |"
	echo "| 1 - Rodar testes                  |"
	echo "| 2 - Limpar os arquivos txt        |"
	echo "| 3 - Limpar os arquivos .gnu       |"
    echo "| 4 - Limpar os arquivos .jpeg      |"
    echo "| 5 - Remover as pastas             |"
	echo "| 6 - Listar arquivos               |"
	# echo "| 6 - Instalar projeto do GitLab    |"
	echo "| 7 - Finalizar Script              |"
	echo " ==================================="

	read opcao

	case $opcao in 
		1)  
			tipoFuncao=(gauss normaL2);
			# tipoFuncao=(eliminacaoGauss gaussJacobi);
            for funcao in "${tipoFuncao[@]}"
            do
                #verifica se o diretório da função já foi criado
                if [ -e "./$funcao" ]
                then
                    echo "Diretório já existe";
                else
                    echo "Criando o diretório";
                    mkdir ./$funcao;
                fi
                #percorre no vetor o que eu quero medir do trabalho 
                dado=L2CACHE;
                valores=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000);
                #coloca em um vetor todos os tamanhos de matriz que eu quero
                echo "Medindo $dado da funcao $funcao, aguarde ...";
                #crio um arquivo com o nome dado e do grupo para buscar as informacoes
                touch dados$dado$funcao.txt
                #crio um arquivo com o nome dado e do grupo para colocar as informações filtradas
                touch plot$dado$funcao.txt
                #compilo com ifdef para não precisar ter + de 1 arquivo.c
               make -B $funcao
                for var in "${valores[@]}" 
                do
                	echo "tamanho = $var na função $funcao testando o dado $dado";    		 	
        	        #roda o likwid e joga os dados num arquivo txt separado 
	                likwid-perfctr -C 3 -g $dado -m -f ./pdeSolver -nx $var -ny $var -i 10 -o stdout$dado.txt >> dados$dado$funcao.txt
                done

                echo "Gerando arquivo plot da função $funcao do dado $dado";
                python3 script.py dados$dado$funcao.txt >> plot$dado$funcao.txt
                touch gnu$dado$funcao.gnu
                echo "set terminal jpeg" >> gnu$dado$funcao.gnu
                echo "set title 'Função $funcao testando $dado'" >> gnu$dado$funcao.gnu
                echo "set key above" >> gnu$dado$funcao.gnu
                echo "set xlabel 'TAMANHO'" >> gnu$dado$funcao.gnu
                echo "set ylabel 'INDICADOR DO TESTE.'" >> gnu$dado$funcao.gnu
                echo "set output 'imagem$dado$funcao.jpeg'" >> gnu$dado$funcao.gnu
                echo "plot 'plot$dado$funcao.txt' title 'data cache miss ratio' with linespoints ls 1" >> gnu$dado$funcao.gnu
                gnuplot -e "load 'gnu$dado$funcao.gnu'; exit"
                mv imagem$dado$funcao.jpeg $funcao/	
                
                dado=L3;
                valores=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000);
                #coloca em um vetor todos os tamanhos de matriz que eu quero
                echo "Medindo $dado da funcao $funcao, aguarde ..."
                #crio um arquivo com o nome dado e do grupo para buscar as informacoes
                touch dados$dado$funcao.txt
                #crio um arquivo com o nome dado e do grupo para colocar as informações filtradas
                touch plot$dado$funcao.txt
                #compilo com ifdef para não precisar ter + de 1 arquivo
                for var in "${valores[@]}" 
                do
                    echo "tamanho = $var na função $funcao testando o dado $dado";  	
                    #roda o likwid e joga os dados num arquivo txt separado 
                    likwid-perfctr -C 3 -g $dado -m -f ./pdeSolver -nx $var -ny $var -i 10 -o stdout$dado.txt>> dados$dado$funcao.txt
                done    

                echo "Gerando arquivo plot da função $funcao do dado $dado";
                python3 script.py dados$dado$funcao.txt >> plot$dado$funcao.txt
                touch gnu$dado$funcao.gnu
                echo "set terminal jpeg" >> gnu$dado$funcao.gnu
                echo "set title 'Função $funcao testando $dado'" >> gnu$dado$funcao.gnu
                echo "set key above" >> gnu$dado$funcao.gnu
                echo "set xlabel 'TAMANHO'" >> gnu$dado$funcao.gnu
                echo "set ylabel 'INDICADOR DO TESTE.'" >> gnu$dado$funcao.gnu
                echo "set output 'imagem$dado$funcao.jpeg'" >> gnu$dado$funcao.gnu
                echo "plot 'plot$dado$funcao.txt' title 'Memory bandwidth [MBytes/s]' with linespoints ls 1" >> gnu$dado$funcao.gnu
                gnuplot -e "load 'gnu$dado$funcao.gnu'; exit"
                mv imagem$dado$funcao.jpeg $funcao/	

                dado=FLOPS_DP;
                valores=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000);
                #coloca em um vetor todos os tamanhos de matriz que eu quero
                echo "Medindo $dado da funcao $funcao, aguarde ..."
                #crio um arquivo com o nome dado e do grupo para buscar as informacoes
                touch dados$dado$funcao.txt
                #crio um arquivo com o nome dado e do grupo para colocar as informações filtradas
                touch plot$dado$funcao.txt
                touch plotAVX$dado$funcao.txt
                #compilo com ifdef para não precisar ter + de 1 arquivo
                for var in "${valores[@]}" 
                do
                    echo "tamanho = $var na função $funcao testando o dado $dado";  	
                    #roda o likwid e joga os dados num arquivo txt separado 
                    likwid-perfctr -C 3 -g $dado -m -f ./pdeSolver -nx $var -ny $var -i 10 -o stdout$dado.txt >> dados$dado$funcao.txt
                done
                python3 script.py dados$dado$funcao.txt DP >> plot$dado$funcao.txt
                echo 'Rodei pro FLOPS_DP'
		      python3 script.py dados$dado$funcao.txt AVX >> plotAVX$dado$funcao.txt
                echo 'Rodei pro AVX FLOPS_DP'

                echo "Gerando arquivo plot da função $funcao do dado $dado";

                touch gnu$dado$funcao.gnu
                echo "set terminal jpeg" >> gnu$dado$funcao.gnu
                echo "set title 'Função $funcao testando $dado'" >> gnu$dado$funcao.gnu
                echo "set key above" >> gnu$dado$funcao.gnu
                echo "set xlabel 'TAMANHO'" >> gnu$dado$funcao.gnu
                echo "set ylabel 'INDICADOR DO TESTE.'" >> gnu$dado$funcao.gnu
                echo "set output 'imagem$dado$funcao.jpeg'" >> gnu$dado$funcao.gnu
                echo "plot 'plot$dado$funcao.txt' title 'FLOPS DP' with linespoints ls 1, 'plotAVX$dado$funcao.txt' title 'AVX FLOPS DP' with linespoints ls 2"  >> gnu$dado$funcao.gnu
                gnuplot -e "load 'gnu$dado$funcao.gnu'; exit"
                mv imagem$dado$funcao.jpeg $funcao/	

                dado=TEMPO;
                valores=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000);
                touch plot$dado$funcao.txt 
                make purge
                make $funcao$dado
                for var in "${valores[@]}" 
                do
                    echo "tamanho = $var na função $funcao testando o dado $dado";  	
                    #roda o likwid e joga os dados num arquivo txt separado 
                    ./pdeSolver -nx $var -ny $var -i 10 -o stdout$dado.txt >> plot$dado$funcao.txt
                done
                echo "Gerando arquivo plot da função $funcao do dado $dado";
                touch gnu$dado$funcao.gnu
                echo "set terminal jpeg" >> gnu$dado$funcao.gnu
                echo "set title 'Função $funcao testando $dado'" >> gnu$dado$funcao.gnu
                echo "set key above" >> gnu$dado$funcao.gnu
                echo "set xlabel 'TAMANHO'" >> gnu$dado$funcao.gnu
                echo "set ylabel 'INDICADOR DO TESTE.'" >> gnu$dado$funcao.gnu
                echo "set output 'imagem$dado$funcao.jpeg'" >> gnu$dado$funcao.gnu
                echo "plot 'plot$dado$funcao.txt' title 'Tempo em milisegundos' with linespoints ls 1" >> gnu$dado$funcao.gnu
                gnuplot -e "load 'gnu$dado$funcao.gnu'; exit"
                mv imagem$dado$funcao.jpeg $funcao/	
            done
            ;;
        2)
            rm *.txt
            ;;
        3)
            rm *.gnu
            ;;
        4)
            rm ptrVet/*.jpeg
            rm rowVet/*.jpeg
            rm colVet/*.jpeg
            ;;
        5)
            rm -r */
            ;;
        6) 
            ls
            ;;

        7) 
            echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
            echo "MODO: POWERSAVE"
            echo "Script Adaptado"
            ;;
    esac
done
