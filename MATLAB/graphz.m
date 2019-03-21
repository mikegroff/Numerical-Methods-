function[retval] = graphz(a,c1,b,c,c3)
hold on;
set(gcf,'Color', 'w');
x = 1:c1; xxx = 1:c3;

plot(xxx,b(xxx),xxx,c(xxx));
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos);
set(ax2,'XAxisLocation','top'); set(ax2,'YAxisLocation','right');
line(x,a(x),'Parent',ax2);

title('Constants of Motion');
xlabel('t'); ylabel('value');
legend('M constant NLS1','M constant NLS2', 'M constant C-N', 'N constant NLS1', 'N constant NLS2','N constant C-N');
hold off;
retval = 0;
end