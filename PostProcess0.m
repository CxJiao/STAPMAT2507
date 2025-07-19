function Post = PostProcess
    Post.Truss = @PostProcess_Truss;
    Post.Plane = @PostProcess_Plane;
    Post.Solid = @PostProcess_Solid;
end

function [] = PostProcess_Truss(Node_X, Node_Y, Dis_X, Dis_Y, beta, Stress, ELE )
    figure();
    hold on;
    for i = 1:size(ELE,1)
        plot([Node_X(ELE(i,1)), Node_X(ELE(i,2))], [Node_Y(ELE(i,1)), Node_Y(ELE(i,2))], 'black', 'LineWidth', 2);
    end
    axis equal;
    
    figure();
    hold on;
    BlockColor = 1000;
    jet_color = colormap(jet(BlockColor));
    [maxValue, ~] = max(Stress);
    [minValue, ~] = min(Stress);

    Node_X = Node_X + beta.*Dis_X;
    Node_Y = Node_Y + beta.*Dis_Y;

    for i = 1:size(ELE,1)
        IndexColor = round((Stress(i)-minValue)*(BlockColor-1)/(maxValue-minValue))+1;
        MappingColor = jet_color(IndexColor,:);
        plot([Node_X(ELE(i,1)), Node_X(ELE(i,2))], [Node_Y(ELE(i,1)), Node_Y(ELE(i,2))], 'color', MappingColor, 'LineWidth', 2);
    end
    axis equal;
    c = colorbar;
    c.Ticks = [0.0 0.2 0.4 0.6 0.8 1.0];
    c.TickLabels = {num2str(minValue), ...
                    num2str((minValue+(maxValue-minValue)/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*2/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*3/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*4/5.0)), ...
                    num2str(maxValue)};
end


function [] = PostProcess_Plane(Node_Coord, Stress, ELE)
    figure();
    hold on;
    for i = 1:size(ELE,1)
        plot([Node_Coord(ELE(i,1),1), Node_Coord(ELE(i,2),1)], ...
             [Node_Coord(ELE(i,1),2), Node_Coord(ELE(i,2),2)], 'black', 'LineWidth', 2);
        plot([Node_Coord(ELE(i,2),1), Node_Coord(ELE(i,3),1)], ...
             [Node_Coord(ELE(i,2),2), Node_Coord(ELE(i,3),2)], 'black', 'LineWidth', 2);
        plot([Node_Coord(ELE(i,1),1), Node_Coord(ELE(i,3),1)], ...
             [Node_Coord(ELE(i,1),2), Node_Coord(ELE(i,3),2)], 'black', 'LineWidth', 2);
    end
    axis equal;
    colorbar;
    Num_element = size(ELE);

    BlockColor = 1000;
    jet_color = colormap(jet(BlockColor));
    [maxValue, ~] = max(Stress);
    [minValue, ~] = min(Stress);

    for i = 1:Num_element
        IndexColor = round((Stress(i)-minValue)*(BlockColor-1)/(maxValue-minValue))+1;
        MappingColor = jet_color(IndexColor,:);
        Node_sum=[ELE(i,1)   Node_Coord(ELE(i,1),1)   Node_Coord(ELE(i,1),2);
                  ELE(i,2)   Node_Coord(ELE(i,2),1)   Node_Coord(ELE(i,2),2);
                  ELE(i,3)   Node_Coord(ELE(i,3),1)   Node_Coord(ELE(i,3),2)];
        fill([Node_sum(1,2) Node_sum(2,2)  Node_sum(3,2)],...
             [Node_sum(1,3) Node_sum(2,3)  Node_sum(3,3)], MappingColor);
    end
    c = colorbar;
    c.Ticks = [0.0 0.2 0.4 0.6 0.8 1.0];
    c.TickLabels = {num2str(minValue), ...
                    num2str((minValue+(maxValue-minValue)/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*2/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*3/5.0)), ...
                    num2str((minValue+(maxValue-minValue)*4/5.0)), ...
                    num2str(maxValue)};
end


function [] = PostProcess_Solid(Node_Coord, Stress, Dis, ELE)
    % Output result for post processing using TecPlot
    global cdata;
    global sdata;
    % Open file
    cdata.IOUT = fopen('.\Data\Solid_Stress_out_1.dat', 'w');
    IOUT = cdata.IOUT;
    [~, NNODE, NUME] = size(Stress);
    [NodeNum, ~] = size(Dis);
    Stress_Von = zeros(NNODE, NUME);
    Stress_Von_Node = zeros(NodeNum, 1);
    for i = 1:NNODE
        for j = 1:NUME
        Stress_Von(i,j) = sqrt(Stress(1,i,j)*Stress(1,i,j)+ ...
                               Stress(2,i,j)*Stress(2,i,j)- ...
                               Stress(1,i,j)*Stress(2,i,j)+ ...
                               3.0*(Stress(4,i,j)*Stress(4,i,j)));
        end
    end

    title = '3D Solid Element ';
    for i = 1:NNODE
        for j = 1:NUME
            Stress_Von_Node(ELE(j,i)) = Stress_Von(i,j);
            disp(ELE(j,i));
            disp(Stress_Von(i,j));
        end
    end

    fprintf(IOUT, 'TITLE = "%s"\n', title);
    fprintf(IOUT, 'VARIABLES = "X", "Y", "Z", "Stress_VonMises"\n');
    fprintf(IOUT, 'ZONE T="Mesh", N=%d, E=%d, F=FEPOINT, ET=BRICK\n', NodeNum, NUME);

    for i = 1:NodeNum
        fprintf(IOUT, '%f   %f   %f   %f\n', Node_Coord(i,1), Node_Coord(i,2), Node_Coord(i,3), Stress_Von_Node(i));
    end
    for i = 1:NUME
        fprintf(IOUT, '%d   %d   %d   %d   %d   %d   %d   %d\n', ELE(i,:)');
    end
    fclose(IOUT);
    disp('Finite element mesh data has been written to mesh.dat');
end

