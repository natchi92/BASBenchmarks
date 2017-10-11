function plotFigures(X1, Y,heading)
%  X:  vector of x data
%  Y:  matrix of y data

% Create figure
figure1 = figure;
n       = ceil(size(Y,1)/2);
k       = size(Y,2);
for i = 1:size(Y,1)
        % Create subplot
        subplot1 = subplot(2,n,i,'Parent',figure1);
        hold(subplot1,'on');

        % Create multiple lines using matrix input to plot
        plot1 = plot(X1,Y(i,:),'Parent',subplot1);

        % Create xlabel
        xlabel('Time (Minutes)');

        % Create title
        title(heading{1,i});

        % Create ylabel
        ylabel('Temperature (^oC)');

        box(subplot1,'on');
        grid(subplot1,'on');
         % Set the remaining axes properties
       % xtick = 1:4:k;
       % set(subplot1,'XTick',xtick','XTickLabel',num2str(1:72)');
        % Create legend
       % legend(subplot1,'show');
end
