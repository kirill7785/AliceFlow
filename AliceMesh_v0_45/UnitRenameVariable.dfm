object FormRenameVar: TFormRenameVar
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Rename project variable'
  ClientHeight = 185
  ClientWidth = 417
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 417
    Height = 185
    Caption = 'Select project variable for rename'
    Color = clMoneyGreen
    ParentBackground = False
    ParentColor = False
    TabOrder = 0
    object LabelRenamecandidate: TLabel
      Left = 16
      Top = 53
      Width = 165
      Height = 13
      Caption = 'List of the names project variables'
    end
    object Label1: TLabel
      Left = 16
      Top = 104
      Width = 342
      Height = 13
      Caption = 
        'Please, enter a new name project variable. First simbol must be ' +
        'equal $'
    end
    object Label2: TLabel
      Left = 24
      Top = 24
      Width = 156
      Height = 13
      Caption = 'Conflict variable name is equal ='
    end
    object LabelConflictName: TLabel
      Left = 200
      Top = 24
      Width = 88
      Height = 13
      Caption = 'LabelConflictName'
    end
    object ComboBox1: TComboBox
      Left = 111
      Top = 72
      Width = 290
      Height = 21
      TabOrder = 0
      Text = 'ComboBox1'
    end
    object EditNewName: TEdit
      Left = 232
      Top = 130
      Width = 169
      Height = 21
      TabOrder = 1
    end
    object Button1: TButton
      Left = 296
      Top = 157
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = Button1Click
    end
  end
end
