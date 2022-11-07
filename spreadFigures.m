% Find all active figures and spread them over the screen

xPad = 50;
yPad = 50;
xBetween = 400;
yBetween = 600;
xBorder = 2400;
yBorder = 1400;
figureHandles = findall(groot, 'Type', 'figure')';
width = 600;
height = 500;
for figureN = figureHandles
    figure(figureN.Number);
    xPos = xPad + mod((figureN.Number-1)*xBetween, xBorder);
    yPos = yBorder - (yPad + floor((figureN.Number-1)*xBetween/xBorder)*yBetween);
    figureN.Position = [xPos, yPos, width, height];
end
