import cv2 as cv

img = cv.imread("TestMap.png")
gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

thresh = cv.threshold(gray, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)[1]

cnts, hierarcy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
cnts = cnts[0] if len(cnts) == 2 else cnts[1]
img_cnt = cv.drawContours(img, cnts, -1, (0,255,0), 3)
cv.imshow('contours', img_cnt)
cv.waitKey()
# for c in cnts:
#     peri = cv.arcLength(c, True)
#     approx = cv.approxPolyDP(c, 0.015 * peri, True)
#     if len(approx) == 4:
#         x,y,w,h = cv.boundingRect(approx)
#         cv.rectangle(img,(x,y),(x+w,y+h),(36,255,12),2)



# cv.imshow('thresh', thresh)
# cv.imshow('image', img)
# cv.waitKey()